import os
import time
import copy
import sys
sys.path.append("..")

import polish.polish_region_implementation as impl
import utils.file_handler as io
import utils.paf_helper as paf
import utils.fasta_handler as fasta

from structures.region import SimpleRegion
import utils.external_tools as tools

def get_region(fork1, fork2, lengthData, q=True):
    
    tid = fork1.qid if q else fork1.rid
    if tid != (fork2.qid if q else fork2.rid): 
        return None
    
    p = [fork1.get_pos_norm(lengthData, q), 
          fork2.get_pos_norm(lengthData, q)]
    return SimpleRegion(tid, min(p[0], p[1]), max(p[0], p[1]))

def path_to_regions(blockPath, q=None):
    
    regions = [] 
    if len(blockPath) > 1:
    
        for i,fork in enumerate(blockPath[:-1]):
            if fork.after_id() == blockPath[i+1].before_id():
                
                if q is not None and not q == fork.is_switch_query():
                    continue                    
                    
                region = SimpleRegion(fork.after_id(), fork.after_pos(), blockPath[i+1].before_pos())
                regions.append(region)
            
    return regions

def get_high_confidence_regions(consensusFa, blockPath, outdir, seqData):
    
    qSeqName = "query"

    #get high confidence regions
    regions = path_to_regions(blockPath, q=True)
    qDict = dict()
    for i,qr in enumerate(regions):
        qDict[qSeqName + "_part" + str(i)] = seqData[qr.chrom][qr.start:qr.end]
    
    qFa = tools.dict2fasta(qDict, outdir + qSeqName)
    qPaf = tools.align_paf_lenient(consensusFa, qFa, outdir + "q_to_consensus")
    qPafRegions = paf.paf_to_regions(qPaf, q=False)
    io.delete_file(qFa)
    return qPafRegions

    
    
    
def polish_contig(tigId, outdir, seqData, lengthData, param):
    cd = os.getcwd()
    os.chdir(outdir)
    
    info = dict()
    startTime = time.time()
    
    seqNames = ["consensus", "hap1", "hap2"]

    region = SimpleRegion(tigId, 0, lengthData[tigId])
    
    #extract reads
    refReads = impl.fetch_ref_reads(param.REF_ALIGNED_READS, region, outdir)
    queryReads = impl.fetch_query_reads(param.QUERY_ALIGNED_READS, region, outdir) 

    #polish sequence
    hybridFa = fasta.write_fasta(outdir + "unpolished.fasta", {tigId : seqData[tigId]})
    refPolishedFa = impl.polish_ref(hybridFa, region, outdir, param, \
                                    refAlignments=param.REF_ALIGNED_READS, outName="ref_polished")
    
    queryPolishedFa = impl.polish_query(refPolishedFa, region, outdir, param, \
                                        queryReads=queryReads, outName="consensus")

    consensusFa = fasta.rename_single_fasta(queryPolishedFa, seqNames[0], toUpper=True)

    #align reads
    refBam = impl.align_ref(consensusFa, refReads, outdir, param, "ref_to_consensus")
    queryBam = impl.align_query(consensusFa, queryReads, outdir, param)

    #get heterozygous variants and phase reads
    highConfVCF = impl.high_conf_hets(consensusFa, refBam, queryBam, outdir, param)
    refHaploBam, queryHaploBam = impl.phase_consensus(consensusFa, highConfVCF, refBam, queryBam, outdir, param)

    hap1Fa, hap2Fa = consensusFa, consensusFa
    
    finalFa1, finalFa2 = outdir + "hap1.fasta", outdir + "hap2.fasta"
    
    if io.file_exists(finalFa1) and os.path.isfile(finalFa2):
        print("Polished haplotypes found, skipping step")
    else:
    
        niter = 1
        realign = False
        while niter > 0:
            hap1Fa = impl.haplotype_polish_ref(1, hap1Fa, refHaploBam, outdir, param, realign)
            hap2Fa = impl.haplotype_polish_ref(2, hap2Fa, refHaploBam, outdir, param, realign)
            realign=True
            
            if param.KEEP_INTERMEDIATE:
                #optional back align step
                hap1Fa = fasta.rename_single_fasta(hap1Fa, seqNames[1])
                hap2Fa = fasta.rename_single_fasta(hap2Fa, seqNames[2])
                tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refpolished")
                tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refpolished")
            
            hap1Fa = impl.haplotype_polish_query(1, hap1Fa, queryHaploBam, outdir, param, realign)
            hap2Fa = impl.haplotype_polish_query(2, hap2Fa, queryHaploBam, outdir, param, realign)
        
            if param.KEEP_INTERMEDIATE:
                #optional back align step
                hap1Fa = fasta.rename_single_fasta(hap1Fa, seqNames[1])
                hap2Fa = fasta.rename_single_fasta(hap2Fa, seqNames[2])
                tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refquerypolished")
                tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refquerypolished")
    
            niter -=1
        
        io.move_file(hap1Fa, finalFa1)
        io.move_file(hap2Fa, finalFa2)

        finalFa1 = fasta.rename_single_fasta(finalFa1, seqNames[1])
        hap2Fa = fasta.rename_single_fasta(hap2Fa, seqNames[2])
        
    
    '''
    
    refBam,queryBam = polisher.phase_reads(consensusFa, highConfVCF, refBam, queryBam, outdir, param)

    #call variants
    refCallsVCF = polisher.call_variants_ref(consensusFa, refBam, outdir, param)
    queryCallsVCF = polisher.call_variants_query(consensusFa, queryBam, outdir, param)

    mappableRegions = get_high_confidence_regions(consensusFa, blockPath, outdir, seqData)
    consensusVariants = polisher.call_variants(consensusFa, queryCallsVCF,refCallsVCF, mappableRegions, outdir)    
    finalVCF = polisher.phase_vcf(consensusFa, consensusVariants, refBam, queryBam, outdir, param)
    '''
    
    fastas = [consensusFa, hap1Fa, hap2Fa]
    graph = tools.construct_graph_msga(fastas, outdir + "graph", normalize=True, renameSeqs=seqNames, baseSeq=seqNames[0])
    tools.index_graph(graph)
    finalVCF = tools.graph_to_vcf(graph, "consensus", "hap", outdir + "graph")
    
    info["time"] = round(time.time() - startTime, 1)
    
    if not param.KEEP_INTERMEDIATE:
        impl.clean_directory(outdir)
    
    #TODO: CLEAN FINAL VCF, VALIDATE HAPLOTYPE, PHASESETS
 
    if param.ANALYSIS:
        impl.analyze_contig(hybridFa, consensusFa, hap1Fa, hap2Fa, tigId, outdir, param)
    
    os.chdir(cd)

    return finalVCF

























def polish_block(blockPath, seqData, lengthData, param):
    cd = os.getcwd()
    os.chdir(param.OUTPUT_DIR)
    
    info = dict()
    startTime = time.time()

    startFork, endFork = blockPath[0], blockPath[-1]
    qRegion = get_region(startFork, endFork, lengthData, q=True)
    rRegion = get_region(startFork, endFork, lengthData, q=False)
    
    seqParts, _ = io.path_to_sequence(blockPath, seqData)
    seq = "".join(seqParts)

    outdir = "./" + "_".join(["blockpolish", str(qRegion.chrom),
                      str(qRegion.start), str(qRegion.end)]) + "/"
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        return

    seqName = "hybrid"
    seqNames = ["consensus", "hap1", "hap2"]

    #extract reads
    refReads = polisher.fetch_ref_reads(param.REF_ALIGNED_READS, rRegion, outdir)
    return
    queryReads = polisher.fetch_query_reads(param.QUERY_ALIGNED_READS, qRegion, outdir) 

    #polish sequence
    hybridFa = tools.dict2fasta({seqName:seq}, outdir + "unpolished")
    refPolishedFa = polisher.polish_ref_reads(hybridFa, refReads, rRegion, seqData, outdir, param)
    consensusFa = polisher.polish_query_reads(refPolishedFa, queryReads, qRegion, seqData, outdir, param)
    consensusFa = tools.rename_fasta(consensusFa, seqNames[0])
    
    #align reads
    refBam = polisher.align_to_consensus_ref(consensusFa, refReads, rRegion, seqData, outdir, param)
    queryBam = polisher.align_to_consensus_query(consensusFa, queryReads, qRegion, seqData, outdir, param)

    #get heterozygous variants and phase reads
    highConfVCF = polisher.high_conf_hets(consensusFa, refBam, queryBam, outdir, param)
    refHaploBam, queryHaploBam = polisher.phase_consensus(consensusFa, highConfVCF, refBam, queryBam, outdir, param)

    hap1Fa, hap2Fa = consensusFa, consensusFa
    
    niter = 1
    realign = False
    while niter > 0:
        hap1Fa = polisher.haplotype_polish_ref(1, hap1Fa, refHaploBam, outdir, param, realign)
        hap2Fa = polisher.haplotype_polish_ref(2, hap2Fa, refHaploBam, outdir, param, realign)
        realign=True
        
        #optional back align step
        hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
        hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])
        tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refpolished")
        tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refpolished")
        
        hap1Fa = polisher.haplotype_polish_query(1, hap1Fa, queryHaploBam, outdir, param, realign)
        hap2Fa = polisher.haplotype_polish_query(2, hap2Fa, queryHaploBam, outdir, param, realign)
    
        #optional back align step
        hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
        hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])
        tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refquerypolished")
        tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refquerypolished")

        niter -=1
    
    hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
    hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])

    
    '''
    
    refBam,queryBam = polisher.phase_reads(consensusFa, highConfVCF, refBam, queryBam, outdir, param)

    #call variants
    refCallsVCF = polisher.call_variants_ref(consensusFa, refBam, outdir, param)
    queryCallsVCF = polisher.call_variants_query(consensusFa, queryBam, outdir, param)

    mappableRegions = get_high_confidence_regions(consensusFa, blockPath, outdir, seqData)
    consensusVariants = polisher.call_variants(consensusFa, queryCallsVCF,refCallsVCF, mappableRegions, outdir)    
    finalVCF = polisher.phase_vcf(consensusFa, consensusVariants, refBam, queryBam, outdir, param)
    '''
    
    fastas = [consensusFa, hap1Fa, hap2Fa]
    graph = tools.construct_graph_msga(fastas, outdir + "graph", normalize=True, renameSeqs=seqNames, baseSeq=seqNames[0])
    tools.index_graph(graph)
    finalVCF = tools.graph_to_vcf(graph, "consensus", "hap", outdir + "graph")
    
    
    info["time"] = round(time.time() - startTime,1)
    
    #polisher.clean_directory(outdir)
    
    #TODO: CLEAN FINAL VCF, VALIDATE HAPLOTYPE, PHASESETS
 
    os.chdir(cd)

    return finalVCF


def polish_interblock(blockPathBefore, blockPathAfter, seqData, lengthData, param):
    
    info = dict()
    startTime = time.time()

    startFork, endFork = blockPathBefore[-1], blockPathAfter[0]
    qRegion = get_region(startFork, endFork, lengthData, q=True)
    rRegion = get_region(startFork, endFork, lengthData, q=False)
    
    seqParts, _ = io.path_to_sequence([startFork, endFork], seqData)
    seq = "".join(seqParts)

    outdir = param.OUTPUT_DIR + "/" + "_".join(["blockpolish", str(qRegion.chrom),
                      str(qRegion.start), str(qRegion.end)]) + "/"
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    seqName = "hybrid"
    seqNames = ["consensus", "hap1", "hap2"]

    #extract reads
    refReads = polisher.fetch_ref_reads(param.REF_ALIGNED_READS, rRegion, outdir)

    #polish sequence x2
    hybridFa = tools.dict2fasta({seqName:seq}, outdir + "unpolished")
    refPolishedFa = polisher.polish_ref_reads(hybridFa, refReads, rRegion, seqData, outdir, param)
    consensusFa = polisher.polish_ref_reads(refPolishedFa, refReads, rRegion, seqData, outdir, param,
                   outName="consensus")
    consensusFa = tools.rename_fasta(consensusFa, seqNames[0])

    #align reads
    refBam = polisher.align_to_consensus_ref(consensusFa, refReads, rRegion, seqData, outdir, param)

    #get heterozygous variants and phase reads
    refHaploBam = polisher.phase_reads_refonly(consensusFa, refBam, outdir, param)
    hap1Fa, hap2Fa = consensusFa, consensusFa
    
    niter = 2
    realign = False
    while niter > 0:
        hap1Fa = polisher.haplotype_polish_ref(1, hap1Fa, refHaploBam, outdir, param, realign)
        hap2Fa = polisher.haplotype_polish_ref(2, hap2Fa, refHaploBam, outdir, param, realign)
        realign=True
        
        #optional back align step
        hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
        hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])
        tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refpolished")
        tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refpolished")
        
        niter -=1
        
    hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
    hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])

    #todo: is phasing good?

    fastas = [consensusFa, hap1Fa, hap2Fa]
    graph = tools.construct_graph_msga(fastas, outdir + "graph", normalize=True, renameSeqs=seqNames, baseSeq=seqNames[0])
    finalVCF = tools.graph_to_vcf(graph, "consensus", "hap", outdir + "graph")
    
    info["time"] = round(time.time() - startTime,1)
    
    #polisher.clean_directory(outdir)
    #TODO: CLEAN FINAL VCF, VALIDATE HAPLOTYPE, PHASESETS

    return finalVCF

def polish_interMblock(blockPathBefore, blockPathAfter, seqData, lengthData, param):
    
    info = dict()
    startTime = time.time()

    startFork, endFork = copy.deepcopy(blockPathBefore[-1]),  copy.deepcopy(blockPathAfter[0])
    startFork.switch_query(); endFork.switch_reference()
    qRegion = get_region(startFork, endFork, lengthData, q=True)
    
    posL = startFork.get_pos_norm(lengthData, q=False)
    ridL = startFork.rid
    rRegionLeft = SimpleRegion(ridL, max(0, posL-1000), min(lengthData[ridL], posL+1000))
   
    posR = endFork.get_pos_norm(lengthData, q=False)
    ridR = endFork.rid
    rRegionRight = SimpleRegion(ridR, max(0, posR-1000), min(lengthData[ridR], posR+1000))

    seqParts, _ = io.path_to_sequence([startFork, endFork], seqData)
    seqN = "".join(seqParts)
    seq = "NNNNNNNNNN".join([x for x in seqN.split("N") if len(x) > 0])

    outdir = param.OUTPUT_DIR + "/" + "_".join(["gappolish", str(qRegion.chrom),
                      str(qRegion.start), str(qRegion.end)]) + "/"
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    seqName = "gap"
    seqNames = ["consensus", "hap1", "hap2"]

    gapFa = tools.dict2fasta({seqName:seq}, outdir + "gap")
    gapBam = tools.align_pacbio(gapFa, param.REF_ALIGNED_READS, outdir + "pacbio_aligned_gap_")

    #extract reads
    refReadsL = polisher.fetch_ref_reads(param.REF_ALIGNED_READS, rRegionLeft, outdir, fileSuffix="_L")
    refReadsR = polisher.fetch_ref_reads(param.REF_ALIGNED_READS, rRegionRight, outdir, fileSuffix="_R")

    niter = 10
    realign = False
    def align(fa):
        refBamL = polisher.align_across_gap_ref(fa, refReadsL, rRegionLeft, seqData, outdir, param, fileSuffix="_L")
        refBamR = polisher.align_across_gap_ref(fa, refReadsR, rRegionRight, seqData, outdir, param,fileSuffix="_R")
        refAlignments = tools.samtools_fetch(refBamL) + tools.samtools_fetch(refBamR)
        refBam = tools.samtools_write(refAlignments, outdir + "ref_aligned_gap", refBamL)
        return refBam
    
    while niter > 0:
        print(niter)
        #align reads
        refBam = align(gapFa)
        consensusFa = polisher.polish_across_gap(gapFa, refBam, seqData, outdir, param, outName="ref_polished")
        consensusFa = tools.rename_fasta(consensusFa, seqNames[0])
        tools.align_pacbio(gapFa, consensusFa, outdir + "consensus_backaligned_gap")
        
        align(consensusFa)

        gapFa=consensusFa
        
        niter -=1


    #align reads
    refBam = polisher.align_to_consensus_ref(consensusFa, refReads, rRegion, seqData, outdir, param)

    #get heterozygous variants and phase reads
    refHaploBam = polisher.phase_reads_refonly(consensusFa, refBam, outdir, param)
    hap1Fa, hap2Fa = consensusFa, consensusFa
    
    niter = 2
    realign = False
    while niter > 0:
        hap1Fa = polisher.haplotype_polish_ref(1, hap1Fa, refHaploBam, outdir, param, realign)
        hap2Fa = polisher.haplotype_polish_ref(2, hap2Fa, refHaploBam, outdir, param, realign)
        realign=True
        
        #optional back align step
        hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
        hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])
        tools.align_pacbio(consensusFa, hap1Fa, outdir + "hap1_backaligned_refpolished")
        tools.align_pacbio(consensusFa, hap2Fa, outdir + "hap2_backaligned_refpolished")
        
        niter -=1
        
    hap1Fa = tools.rename_fasta(hap1Fa, seqNames[1])
    hap2Fa = tools.rename_fasta(hap2Fa, seqNames[2])

    #todo: is phasing good?

    fastas = [consensusFa, hap1Fa, hap2Fa]
    graph = tools.construct_graph_msga(fastas, outdir + "graph", normalize=True, renameSeqs=seqNames, baseSeq=seqNames[0])
    finalVCF = tools.graph_to_vcf(graph, "consensus", "hap", outdir + "graph")
    
    info["time"] = round(time.time() - startTime,1)
    
    #polisher.clean_directory(outdir)
    #TODO: CLEAN FINAL VCF, VALIDATE HAPLOTYPE, PHASESETS

    return finalVCF
