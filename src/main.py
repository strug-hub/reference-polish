import random, math, re
import parameters
import utils.file_handler as io
import utils.fasta_handler as fasta
import utils.external_tools as tools

import structures.region as rgn
import ref_polish_helper as helper
#import analysis.ref_polish_analysis as analyzer
#import variants as var

'''
def make_haplo_vcfs():
    
    param = parameters.get_parameters_reference_polish()

    print("Parsing parameters...")
    seqData = dict()

    print("Reading reference fasta...")
    refData = fasta.read_fasta(param.REF_FA)
    seqData.update(refData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}

    #isolate region of interest
    region = rgn.region_from_string(param.REF_REGION, lengthData)


    for CFID in [ \
            "CF001", "CF002", "CF003", "CF004", "CF006", "CF007", "CF010", "CF011", "CF013", "CF014", "CF016", 
            "CF022", "CF024", "CF045", "CF047", "CF049", "CF052", "CF060", "CF062", "CF063", "CF066", "CF067",
            "CF071", "CF072", "CF073", "CF075", "CF076", "CF077", "CF078"]:

        param = parameters.get_parameters_reference_polish(CFID)  
        prefix = helper.file_prefix(region, param)
        regionFa, alignedBam = helper.isolate_region(seqData, param.REF_ALIGNED_READS, region, prefix)
        hapAFa = prefix + "hapA_final.polished.fasta"
        hapBFa = prefix + "hapB_final.polished.fasta"
            
        pafFileA = tools.align_paf_very_lenient(regionFa, hapAFa, prefix + "_hapA")
        pafFileB = tools.align_paf_very_lenient(regionFa, hapBFa, prefix + "_hapB")

        variantSetA = paf_helper.get_variantset(pafFileA, regionFa, hapAFa)
        variantSetB = paf_helper.get_variantset(pafFileB, regionFa, hapBFa)
        
        variantSet = var.combine_as_genotypes(variantSetA, variantSetB, phased=True)
        vcf = variantSet.write_vcf(prefix + "_haps.vcf")
        variantSet.shift_all(region.start, region.chrom)
        vcfRef = variantSet.write_vcf(prefix + "_hg38_haps.vcf")


        phasedVcf = tools.whatshap_phase([alignedBam], vcf, regionFa, prefix, indels=False)
        variantSet = var.vcf_to_variantcallset(phasedVcf)

        variantSet.shift_all(region.start, region.chrom)
        variantSet.write_vcf(phasedVcf)

        #vcfA = variantSetA.write_vcf(prefix + "_hapA.vcf")
        #vcfB = variantSetB.write_vcf(prefix + "_hapB.vcf")

        #vcfAPhased = tools.whatshap_phase([alignedBam], vcfA, regionFa, prefix + "_A", indels=False)
        #vcfBPhased = tools.whatshap_phase([alignedBam], vcfB, regionFa, prefix + "_B", indels=False)
        
        #whVariantSetA = var.vcf_to_variantcallset(vcfAPhased)
        #whVariantSetB = var.vcf_to_variantcallset(vcfBPhased)

        #whVariantSetA.shift_all(region.start, region.chrom)
        #whVariantSetB.shift_all(region.start, region.chrom)

        #vcfA = whVariantSetA.write_vcf(prefix + "_A.whatshap.phase.vcf")
        #vcfB = whVariantSetB.write_vcf(prefix + "_B.whatshap.phase.vcf")
''' 
    
def main_human_reference_polish():
    
    param = parameters.get_parameters()
    
    print("Parsing parameters...")
    seqData = dict()

    print("Reading reference fasta...")
    refData = fasta.read_fasta(param.REF_FA)
    seqData.update(refData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}

    backAlign = True
    plot = True

    #isolate region of interest
    targetRegion = rgn.region_from_string(param.REF_REGION, lengthData)

    '''
    for CFID in [ \
            #"CF001", "CF002", "CF003",
             "CF004", "CF006", "CF007", "CF010", "CF011", "CF013", "CF014", "CF016", 
            "CF022", "CF024", "CF045", "CF047", "CF049", "CF052", "CF060", "CF062", "CF063", "CF066", "CF067",
            "CF071", "CF072", "CF073", "CF075", "CF076", "CF077", "CF078"]:
                
        param = parameters.get_parameters_reference_polish(CFID)
    ''' 
    
    faDict, faDictOrder, bamDict = dict(), [], dict()

    print(param.ID)
    
    outdir = param.OUTPUT_DIR
    
    if param.DOWNSAMPLE is not None and param.DOWNSAMPLE < 1:
        outdir = outdir + "_downsampled_" + str(param.DOWNSAMPLE)
    
    io.make_dir(outdir)
    chromName = re.sub('[^a-zA-Z0-9 \n\.]', '-', targetRegion.chrom)
    prefix = outdir + "/" + "_".join([chromName, str(targetRegion.start), str(targetRegion.end)])

    # get sequence
    seqDict = dict()
    fid = "_".join([targetRegion.chrom, str(targetRegion.start), str(targetRegion.end)])
    seqDict[fid] = seqData[targetRegion.chrom][targetRegion.start:targetRegion.end]
    regionFa = fasta.write_fasta(prefix, seqDict, index=True)
        
    faDict["HG38"] = regionFa
    faDictOrder.append("HG38")

    bamFile = param.READS

    # downsample
    if param.DOWNSAMPLE is not None and param.DOWNSAMPLE < 1:
        print("hello")
        if io.file_exists(prefix + "_downsampled.bam"):
            bamFile = prefix + "_downsampled.bam"
        else:
            alignDict = tools.samtools_fetch(param.READS, asDict=True)
            qnames = [x for x in alignDict]
            random.shuffle(qnames)
            half = math.ceil(len(qnames) * param.DOWNSAMPLE)
            keepNames = set(qnames[:half])
    
            filteredAlns = []
            for qname in keepNames: filteredAlns.extend(alignDict[qname])
            bamFile = tools.samtools_write(filteredAlns, prefix + "_downsampled", param.READS)


    # polish with all reads
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

    faDict["Polished"] = polishedFa
    faDictOrder.append("Polished")

    # phase reads
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
        
    
    phaseblocks = helper.phaseblock_split(phasedVCF)
    
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    
    bamDict["HapA"] = bamA
    bamDict["HapB"] = bamB
    bamDict["Unphased"] = bamUnphased

    # haplotype polish

    def haplo_polish(fa, hapBam, prefix, niter=1, chunkSize=None):
        
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

    if io.file_exists(prefix + "hapA.polished.fasta"):
        faA = prefix + "hapA.polished.fasta"
    else:
        faA_ = haplo_polish(polishedFa, bamA, prefix + "hapA_inter", niter=2, chunkSize=5000)
        faA = haplo_polish(faA_, bamA, prefix + "hapA", niter=2)

    faDict["Haplotype A"] = faA
    faDictOrder.append("Haplotype A")

    if io.file_exists(prefix + "hapB.polished.fasta"):
        faB = prefix + "hapB.polished.fasta"
    else:
        faB_ = haplo_polish(polishedFa, bamB, prefix + "hapB_inter", niter=2, chunkSize=5000)
        faB = haplo_polish(faB_, bamB, prefix + "hapB", niter=2)

    faDict["Haplotype B"] = faB
    faDictOrder.append("Haplotype B")

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
        
    def clean_vcf(vcfFile, pblocks, vcfOut):
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
    
    vcfA = clean_vcf(vcfA_, refBlocksA, prefix + "_hap1.vcf")
    vcfB = clean_vcf(vcfB_, refBlocksB, prefix + "_hap2.vcf")

    
    #if False:
    tools.align_pacbio(param.REF_FA, bamA, prefix + "_h1_hg38")
    tools.align_pacbio(param.REF_FA, bamB, prefix + "_h2_hg38")
    tools.align_pacbio(param.REF_FA, bamUnphased, prefix + "_unphased_hg38")

    tools.align_pacbio(backFa, polishedFa, prefix + "_backaligned_polishAll")
    tools.align_pacbio(backFa, polishedFaA, prefix + "_backaligned_polishA1")
    tools.align_pacbio(backFa, polishedFaA2, prefix + "_backaligned_polishA2")
    tools.align_pacbio(backFa, polishedFaB, prefix + "_backaligned_polishB1")
    tools.align_pacbio(backFa, polishedFaB2, prefix + "_backaligned_polishB2")


'''
    if backAlign:
        backFa = param.REF_FA
        #backFa = regionFa
        tools.align_pacbio(backFa, polishedFa, prefix + "_backaligned_polishAll")
        tools.align_pacbio(backFa, polishedFaA, prefix + "_backaligned_polishA1")
        tools.align_pacbio(backFa, polishedFaA2, prefix + "_backaligned_polishA2")
        tools.align_pacbio(backFa, polishedFaB, prefix + "_backaligned_polishB1")
        tools.align_pacbio(backFa, polishedFaB2, prefix + "_backaligned_polishB2")


if plot:
    analyzer.plot_alignments(faDict, bamDict, prefix, faDictOrder=faDictOrder, image="png")
'''
    
        
    

if __name__== "__main__":
  #make_haplo_vcfs() 
  main_human_reference_polish()
  print("done")
  #exit()