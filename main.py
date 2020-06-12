import random, math, re
import parameters
import vcf
import utils.file_handler as io
import utils.fasta_handler as fasta
import utils.external_tools as tools
import utils.region as rgn

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
    

def main():
    
    param = parameters.get_parameters()
    print(param.ID)  

    #isolate region of interest
    #------------------------------
    targetRegion = rgn.region_from_string(param.REF_REGION)
    fid = "_".join([targetRegion.chrom, str(targetRegion.start), str(targetRegion.end)])

    targetDict = {}
    targetDict[fid] = fasta.faidx_fetch(param.REF_FA, targetRegion)

    
    outdir = param.OUTPUT_DIR
    if param.DOWNSAMPLE is not None and param.DOWNSAMPLE < 1:
        outdir = outdir + "_downsampled_" + str(param.DOWNSAMPLE)
    io.make_dir(outdir)

    chromName = re.sub('[^a-zA-Z0-9 \n\.]', '-', targetRegion.chrom)
    prefix = outdir + "/" + "_".join([chromName, str(targetRegion.start), str(targetRegion.end)])

    regionFa = fasta.write_fasta(prefix, targetDict, index=True)

    # get reads
    #------------------------------
    
    bamFile = param.READS

    # downsample
    if param.DOWNSAMPLE is not None and param.DOWNSAMPLE < 1:
        bamFile = downsample(param.READS, param.DOWNSAMPLE, prefix + "_downsampled")

    # polish with all reads
    #------------------------------

    niter = 2
    fa = regionFa
    
    if io.file_exists(prefix + ".polished.fasta"):
        polishedFa = prefix + ".polished.fasta"
    else:
        for i in range(1,niter+1):
            
            currentPrefix = prefix + "_iter_" + str(i)
    
            bam = tools.align_pacbio(fa, bamFile, currentPrefix)            
            fa = tools.pacbio_polish(bam, fa, currentPrefix, outputFasta=True, chunkSize=5000)
            vcf = currentPrefix + ".consensus.vcf"
    
            io.delete_file(vcf)
            io.delete_file(bam)
    
        polishedFa = io.rename_file(fa, prefix + ".polished.fasta")
        io.rename_file(fa + ".fai", prefix + ".polished.fasta.fai")

    # phase reads
    #------------------------------
        
    if io.file_exists(prefix + ".longshot.vcf"):
        phasedVCF = prefix + ".longshot.vcf"
    else:

        polishedBam = tools.align_pacbio(polishedFa, param.READS, prefix + "_consensus")
        
        alignDict = tools.samtools_fetch(polishedBam, asDict=True)
        lenDict = dict()
        lenDict = { qname : sum([x.qlen for x in alignDict[qname]]) for qname in alignDict}
    
        LENGTH_THRESH=3000
        
        filteredAlns = []
        for qname in lenDict:
            if lenDict[qname] > LENGTH_THRESH: filteredAlns.extend(alignDict[qname])
        
        filteredBam = tools.samtools_write(filteredAlns, prefix + "_filtered", polishedBam)
        phasedVCF = tools.longshot_genotype(filteredBam, polishedFa, prefix, writeBams=True)
        
    phaseblocks = phaseblock_split(phasedVCF)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)

    # haplotype polish
    #------------------------------
    
    if io.file_exists(prefix + "hapA.polished.fasta"):
        faA = prefix + "hapA.polished.fasta"
    else:
        faA_ = haplo_polish(polishedFa, bamA, bamUnphased, prefix + "hapA_inter", niter=2, chunkSize=5000)
        faA = haplo_polish(faA_, bamA, bamUnphased, prefix + "hapA", niter=2)

    if io.file_exists(prefix + "hapB.polished.fasta"):
        faB = prefix + "hapB.polished.fasta"
    else:
        faB_ = haplo_polish(polishedFa, bamB, bamUnphased, prefix + "hapB_inter", niter=2, chunkSize=5000)
        faB = haplo_polish(faB_, bamB, bamUnphased, prefix + "hapB", niter=2)

    # call variation
    #------------------------------
    
    pafA_ = tools.align_paf_lenient(faA, polishedFa, prefix + "_hap1_topolished")
    pafB_ = tools.align_paf_lenient(faB, polishedFa, prefix + "_hap2_topolished")
    blocksA = tools.paf_liftover(pafA_, phaseblocks, minAlnLen=2500)
    blocksB = tools.paf_liftover(pafB_, phaseblocks, minAlnLen=2500)

    pafA = tools.align_paf_lenient(regionFa, faA, prefix + "_hap1")
    pafB = tools.align_paf_lenient(regionFa, faB, prefix + "_hap2")
    refBlocksA = tools.paf_liftover(pafA, blocksA, minAlnLen=2500)
    refBlocksB = tools.paf_liftover(pafB, blocksB, minAlnLen=2500)

    vcfA_ = tools.paf_call(pafA, regionFa, prefix + "_hap1_temp", vcfSample=param.ID)
    vcfB_ = tools.paf_call(pafB, regionFa, prefix + "_hap2_temp", vcfSample=param.ID)

    #fid, flen = helper.get_fasta_id(polishedFa), helper.get_fasta_len(polishedFa)
    #lengthData[fid] = flen
        
    vcfA = clean_vcf(vcfA_, targetRegion, refBlocksA, prefix + "_hap1.vcf")
    vcfB = clean_vcf(vcfB_, targetRegion, refBlocksB, prefix + "_hap2.vcf")

    print(vcfA)
    print(vcfB)

    '''
    tools.align_pacbio(param.REF_FA, bamA, prefix + "_h1_hg38")
    tools.align_pacbio(param.REF_FA, bamB, prefix + "_h2_hg38")
    tools.align_pacbio(param.REF_FA, bamUnphased, prefix + "_unphased_hg38")

    tools.align_pacbio(backFa, polishedFa, prefix + "_backaligned_polishAll")
    tools.align_pacbio(backFa, polishedFaA, prefix + "_backaligned_polishA1")
    tools.align_pacbio(backFa, polishedFaA2, prefix + "_backaligned_polishA2")
    tools.align_pacbio(backFa, polishedFaB, prefix + "_backaligned_polishB1")
    tools.align_pacbio(backFa, polishedFaB2, prefix + "_backaligned_polishB2")
    '''

if __name__== "__main__":
  main()
