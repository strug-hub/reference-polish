import random, math
import vcf
import file_handler as io
import external_tools as tools
import region as rgn

def downsample(bam, proportion, outPrefix):

    random.seed(1337)
    
    alignDict = tools.samtools_fetch(bam, asDict=True)
    qnames = [x for x in alignDict]
    random.shuffle(qnames)
    split = math.ceil(len(qnames) * proportion)
    keepNames = set(qnames[:split])
    
    filteredAlns = []
    for qname in keepNames: filteredAlns.extend(alignDict[qname])
    bamFile = tools.samtools_write(filteredAlns, outPrefix, bam)
    return bamFile

def phaseblock_split(phasedVCF, chrom=None):
    # use chrom=None if there is only one chrom, else specify chrom name to 
    # only consider variants on that chromosome
    #todo: overlapping blocks??
    
    phaseblocks = dict()
    for v in vcf.Reader(open(phasedVCF, "r")):
        if chrom is not None and v.CHROM != chrom:
            continue
        try:
            ps = v.samples[0]["PS"]
        except:
            ps = None
        if ps not in phaseblocks: phaseblocks[ps] = []
        phaseblocks[ps].append(v)
    
    blocks = []
    for ps in phaseblocks:
        if ps is None or ps == "." : continue
        positions = [v.POS for v in phaseblocks[ps]]
        
        if len(positions) < 2: continue
    
        start, end = min(positions), max(positions)
        blocks.append(rgn.SimpleRegion(phaseblocks[ps][0].CHROM, start, end))
    
    return blocks

def haplo_polish(fa, hapBam, bamUnphased, prefix, niter=1, chunkSize=None):
        
        for i in range(1,niter+1):
        
            currentPrefix = prefix + "_iter_" + str(i)

            reads = tools.samtools_fetch(hapBam)
            unphased = tools.samtools_fetch(bamUnphased)
            random.shuffle(unphased)
            half = math.ceil(len(unphased)/2)
            reads = reads + unphased[:half]
            tempBam = tools.samtools_write(reads, currentPrefix + "_TEMP_assigned", hapBam)
            bam = tools.align_pacbio(fa, tempBam, currentPrefix)
            io.delete_file(tempBam)

            fa = tools.pacbio_polish(bam, fa, currentPrefix, outputFasta=True, chunkSize=chunkSize)
            vcf = currentPrefix + ".consensus.vcf"

            io.delete_file(vcf)
            io.delete_file(bam)
            
        finalFa = io.rename_file(fa, prefix + ".polished.fasta")
        io.rename_file(fa + ".fai", prefix + ".polished.fasta.fai")

        return finalFa
    
def clean_vcf(vcfFile, targetRegion, pblocks, vcfOut):
    writer = open(vcfOut, "w")

    for line in open(vcfFile, "r"):
        if line[0] == "#" or len(line) < 1: 
            writer.writelines(line)
            continue
        
        s = line.rstrip().split("\t")
        
        for region in pblocks:
            if region.contains_pos(int(s[1])):
                s[8] = s[8] + ":PS"
                s[9] = s[9] + ":" + str(region.start)
                break
            
        s[0] = targetRegion.chrom
        s[1] = str(int(s[1]) + targetRegion.start)

        writer.writelines("\t".join(s) + "\n")
    return vcfOut
    