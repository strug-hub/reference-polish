import re
import parameters
import file_handler as io
import external_tools as tools
import region as rgn
import helper

def main():
    
    param = parameters.get_parameters()
    print(param.ID)  

    #isolate region of interest
    #------------------------------
    targetRegion = rgn.region_from_string(param.REF_REGION)
    fid = "_".join([targetRegion.chrom, str(targetRegion.start), str(targetRegion.end)])

    targetDict = {}
    targetDict[fid] = io.faidx_fetch(param.REF_FA, targetRegion)

    
    outdir = param.OUTPUT_DIR
    if param.DOWNSAMPLE is not None and param.DOWNSAMPLE < 1:
        outdir = outdir + "_downsampled_" + str(param.DOWNSAMPLE)
    io.make_dir(outdir)

    chromName = re.sub('[^a-zA-Z0-9 \n\.]', '-', targetRegion.chrom)
    prefix = outdir + "/" + "_".join([chromName, str(targetRegion.start), str(targetRegion.end)])

    regionFa = io.write_fasta(prefix, targetDict, index=True)

    # get reads
    #------------------------------
    
    bamFile = param.READS

    # downsample
    if param.DOWNSAMPLE is not None and param.DOWNSAMPLE < 1:
        bamFile = helper.downsample(param.READS, param.DOWNSAMPLE, prefix + "_downsampled")

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
        
    phaseblocks = helper.phaseblock_split(phasedVCF)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)

    # haplotype polish
    #------------------------------
    
    if io.file_exists(prefix + "hapA.polished.fasta"):
        faA = prefix + "hapA.polished.fasta"
    else:
        faA_ = helper.haplo_polish(polishedFa, bamA, bamUnphased, prefix + "hapA_inter", niter=2, chunkSize=5000)
        faA = helper.haplo_polish(faA_, bamA, bamUnphased, prefix + "hapA", niter=2)

    if io.file_exists(prefix + "hapB.polished.fasta"):
        faB = prefix + "hapB.polished.fasta"
    else:
        faB_ = helper.haplo_polish(polishedFa, bamB, bamUnphased, prefix + "hapB_inter", niter=2, chunkSize=5000)
        faB = helper.haplo_polish(faB_, bamB, bamUnphased, prefix + "hapB", niter=2)

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
        
    vcfA = helper.clean_vcf(vcfA_, targetRegion, refBlocksA, prefix + "_hap1.vcf")
    vcfB = helper.clean_vcf(vcfB_, targetRegion, refBlocksB, prefix + "_hap2.vcf")

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
