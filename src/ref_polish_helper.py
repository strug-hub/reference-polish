import re
import random, math
import vcf

import utils.external_tools as tools
import utils.file_handler as io
import utils.fasta_handler as fasta
from structures.region import SimpleRegion

'''
import helper
import external_tools as tools
import paf_helper as paf
import regions as rgn
import variants as var
'''

def file_prefix(region, param):
    chromName = re.sub('[^a-zA-Z0-9 \n\.]', '-', region.chrom)
    
    
    return param.OUTPUT_DIR + "/" + "_".join([chromName, str(region.start), str(region.end)])

def isolate_region(seqData, bamFile, region, prefix, unaligned=False):
    
    seqDict = dict()
    fid = "_".join([region.chrom, str(region.start), str(region.end)])
    seqDict[fid] = seqData[region.chrom][region.start:region.end]

    isolationFa = fasta.write_fasta(prefix, seqDict, index=True)
    subsetBam = tools.samtools_subset(bamFile, region, prefix + "_TEMP_")
    
    if unaligned:
        unalignedBam = tools.unalign_bam(subsetBam, prefix)
        io.delete_file(subsetBam)
        return [isolationFa, unalignedBam]
    
    realignedBam = tools.align_pacbio(isolationFa, subsetBam, prefix)
    io.delete_file(subsetBam)
    return [isolationFa, realignedBam]

def iterative_polish(fastaFile, bamFile, prefix, niter=1, 
                     bamUnphased=None, backAlign=False, 
                     useMapFilter=False, chunkSize=None, region=None):
         
    i = 1
    currentFastaFile = fastaFile
    currentBamFile = bamFile
    currentUnphasedBamFile = bamUnphased

    while True:
        
        currentPrefix = prefix + "_iter_" + str(i)
        
        assignedBam = currentBamFile
        if bamUnphased is not None:
            reads = tools.samtools_fetch(currentBamFile)
            unphased = tools.samtools_fetch(currentUnphasedBamFile)
            random.shuffle(unphased)
            half = math.ceil(len(unphased)/2)
            reads = reads + unphased[:half]
            assignedBamTemp = tools.samtools_write(reads, currentPrefix + "_TEMP_assigned", bamFile)
            assignedBam = tools.align_pacbio(currentFastaFile, assignedBamTemp, currentPrefix + "_TEMP_aligned")
            io.delete_file(assignedBamTemp)
            
        prevFastaFile = currentFastaFile
        currentFastaFile = tools.pacbio_polish(assignedBam, currentFastaFile, currentPrefix, region=region,
                                               outputFasta=True, useMapFilter=useMapFilter, chunkSize=chunkSize)
        if assignedBam != bamFile:
            io.delete_file(assignedBam)
        if prevFastaFile != fastaFile:
            io.delete_file(prevFastaFile)
        
        # align new fasta back to original to visualize differences
        if backAlign: 
            tools.align_pacbio(fastaFile, currentFastaFile, currentPrefix + "_backaligned")
            
        consensusVCF = currentPrefix + ".consensus.vcf"
        io.delete_file(consensusVCF)


        if i >= niter:
            io.delete_file(consensusVCF)
            break

        prevBamFile = currentBamFile
        currentBamFile = tools.align_pacbio(currentFastaFile, bamFile, currentPrefix)
        if prevBamFile != bamFile:
            io.delete_file(prevBamFile)

        if bamUnphased is not None:
            prevBamFile = currentUnphasedBamFile
            currentUnphasedBamFile = tools.align_pacbio(currentFastaFile, bamUnphased, currentPrefix + "_unphased")
            if prevBamFile != bamUnphased:
                io.delete_file(prevBamFile)

        i +=1

    if currentBamFile != bamFile:
        io.delete_file(currentBamFile)
    if bamUnphased is not None and currentUnphasedBamFile != bamUnphased:
        io.delete_file(currentUnphasedBamFile)
        
    finalPolishedFa = io.rename_file(currentFastaFile, prefix + ".polished.fasta")
    io.rename_file(currentFastaFile + ".fai", prefix + ".polished.fasta.fai")

    return finalPolishedFa

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
        blocks.append(SimpleRegion(phaseblocks[ps][0].CHROM, start, end))
    
    return blocks


'''

def polish_structural_haplotigs(mainContig, svContig, seqData, lengthData, param, dummy=False):
        
    tigFa = helper.contig_fasta(mainContig, seqData, param)
    svTigFa = helper.contig_fasta(svContig, seqData, param)

    (s1,e1), (s2,e2) = paf.get_alignment_bounds(tigFa, svTigFa, param.OUTPUT_DIR + "/" + svContig)
    svRegion = rgn.SimpleRegion(svContig, 0, lengthData[svContig]-1)
    svBuffer = (s2, len(svRegion) - e2)
    mainRegion = rgn.SimpleRegion(mainContig, s1 - svBuffer[0], e1 + svBuffer[1])

    segment = rgn.GraphSegment(mainRegion, svRegion)

    # Isolate and get reads
    prefix = helper.file_prefix(mainRegion, param)
    regionFa, regionBam = helper.isolate_region(mainRegion, seqData, param)  
    svPrefix = helper.file_prefix(svRegion, param)
    svReads = tools.samtools_subset(param.REF_ALIGNED_READS, svRegion, svPrefix)

    if dummy:
        segment.set_sequence(regionFa, svTigFa)
        print("Resolving structural haplotigs:", mainContig, svContig)
        return segment
    
    # Phase all reads
    mainBam, svBam, bamUnphased = helper.phase_sv_reads(regionFa, regionBam, svTigFa, svReads, \
                                                    param.OUTPUT_DIR + "/resolve_structural_haplotigs", param)

    segment.set_reads(mainBam, svBam, bamUnphased)

    bamA = tools.align_pacbio(mainBam, regionFa, prefix + "_phasedA")
    bamUnphasedA = tools.align_pacbio(mainBam, regionFa, prefix + "_unphasedA")
    mainPolishedFa = iterative_polish(regionFa, bamA, prefix + "_main", bamUnphased=bamUnphasedA, niter=2)
    
    bamB = tools.align_pacbio(svBam, regionFa, prefix + "_phasedB")
    bamUnphasedB = tools.align_pacbio(svBam, regionFa, prefix + "_unphasedB")
    svPolishedFa = iterative_polish(svTigFa, bamB, prefix + "_sv", bamUnphased=bamUnphasedB, niter=2)

    segment.set_sequence(mainPolishedFa, svPolishedFa)
    
    return segment

def phase_region(fa, reads, region, param, entireRegion=False):
    prefix = helper.file_prefix(region, param)
    regionToUse = region if not entireRegion else None
    longshotVCF = tools.longshot_genotype(reads, fa, prefix, regionToUse, writeBams=True)
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    return (longshotVCF, bamA, bamB, bamUnphased)

def phaseblock_split(phasedVCF, region, lengthData, param, collapseBuffer=300):

    variantSetAll = var.vcf_to_variantcallset(phasedVCF, region)
    variantPhaseSets = var.split_by_ps(variantSetAll)
    
    phasedSegments = []
    for ps in variantPhaseSets:
        if ps == ".": continue
        positions = variantPhaseSets[ps].get_positions()
        isHeterozygous = [variantPhaseSets[ps].get(pos).is_heterozygous() for pos in positions]
        positions = [pos for pos,het in zip(positions, isHeterozygous) if het]

        if len(positions) < 2: continue
    
        start, end = min(positions), max(positions)
        phaseBlock = rgn.SimpleRegion(region.chrom, start, end, lengthData)
        phasedSegment = rgn.GraphSegment(phaseBlock, phaseBlock)
        phasedSegments.append(phasedSegment)
    
    unphasedBlocks = rgn.region_difference(region, phasedSegments)
    unphasedSegments = [rgn.GraphSegment(phaseBlock, None) for phaseBlock in unphasedBlocks]
        
    phasedSegments = sorted(phasedSegments, key=lambda r: r.start)
    unphasedSegments = sorted(unphasedSegments, key=lambda r: r.start)
    if len(unphasedSegments) > 0 and len(phasedSegments) > 0 and \
        unphasedSegments[0].start < phasedSegments[0].start and \
        len(unphasedSegments[0]) < collapseBuffer:
            
            phasedSegments[0].start = region.start
            unphasedSegments.pop(0)
        
    phasedSegments = sorted(phasedSegments, key=lambda r: r.end)
    unphasedSegments = sorted(unphasedSegments, key=lambda r: r.end)
    if len(unphasedSegments) > 0 and len(phasedSegments) > 0 and \
        unphasedSegments[-1].end > phasedSegments[-1].end and \
        len(unphasedSegments[-1]) < collapseBuffer:
            phasedSegments[-1].end = region.end
            unphasedSegments.pop()

    #phasedSegments = [segment.extend(buffer) for segment in phasedSegments]
    #unphasedSegments = [segment.extend(buffer) for segment in unphasedSegments]        

    return (phasedSegments, unphasedSegments)


def polish_region_haplotypes(phasedRegion, bamA, bamB, bamUnphased, seqData, lengthData, param,
                             backAlign=False, niter=1, chunkSize=None, dummy=False):
    
    # Isolate region
    regionFa, regionBamA, regionBamB, regionBamU = \
            helper.isolate_phased_region(phasedRegion, seqData, bamA, bamB, bamUnphased, param)
    prefix = helper.file_prefix(phasedRegion, param)
   
    phasedSegment = rgn.GraphSegment(phasedRegion, phasedRegion)
    phasedSegment.set_reads(regionBamA, regionBamB, regionBamU)

    if dummy:
        phasedSegment.set_sequence(regionFa, regionFa)
        print("Resolving haplotigs:", phasedRegion)
        return phasedSegment
        
    # Polish using all reads first??
    #polishedFa = tools.pacbio_polish(regionBam, regionFa, prefix, outputFasta=True)

    # Polish haplotypes 
    consensusFaA = iterative_polish(regionFa, regionBamA, prefix + "_A", bamUnphased=regionBamU, \
                                    niter=niter, useMapFilter=True, backAlign=backAlign, chunkSize=chunkSize)
    consensusFaB = iterative_polish(regionFa, regionBamB, prefix + "_B", bamUnphased=regionBamU, \
                                    niter=niter, useMapFilter=True, backAlign=backAlign, chunkSize=chunkSize)
    
    #pafBtoA = tools.align_paf_lenient(consensusFaA, consensusFaB, prefix + "_BtoA")
    #variantSetB = paf.get_variantset(pafBtoA, consensusFaA)
        
    #readIdsA = list(set(helper.fetch_read_ids(regionBamA)))
    #readIdsB = list(set(helper.fetch_read_ids(regionBamB)))
    #readIdsU = list(set(helper.fetch_read_ids(regionBamU)))
    
    phasedSegment.set_sequence(consensusFaA, consensusFaB)
    #phasedSegment.add_variants(var.VariantCallSet(), variantSetB = variantSetB)
    
    return phasedSegment


def polish_region_unphased(unphasedRegion, bamA, bamB, bamUnphased, seqData, lengthData, param, dummy=False):

    # Isolate region
    regionFa, regionBamA, regionBamB, regionBamU = \
            helper.isolate_phased_region(unphasedRegion, seqData, bamA, bamB, bamUnphased, param)
    prefix = helper.file_prefix(unphasedRegion, param)
    
    phasedSegment = rgn.GraphSegment(unphasedRegion)
    phasedSegment.set_reads(bamU=regionBamU)

    if dummy:
        phasedSegment.set_sequence(regionFa)
        print("Resolving unphased contig:", unphasedRegion)
        return phasedSegment
    
    # Unphase reads
    bamAll = helper.combine_reads([regionBamA, regionBamB, regionBamU], prefix)
   
    # Polish using all reads
    consensusFa = iterative_polish(regionFa, bamAll, prefix + "_all",  useMapFilter=True, backAlign=False)
    #readIds = list(set(helper.fetch_read_ids(bamAll)))

    phasedSegment.add_fasta(consensusFa, None)
    #phasedSegment.add_variants(var.VariantCallSet(), var.VariantCallSet())

    return phasedSegment


def inject_pbsv_calls(fa, bamFile, prefix, medianFilter=True):
    
    def inject_sv(sequence, callset):
        for pos in reversed(callset.get_positions()):
            variant = callset.get(pos)
            allele = variant.alleleB if variant.alleleA is None else variant.alleleA
            ref, alt = allele.get_ref(), allele.get_alt()
            sequence = sequence[:pos] + alt + sequence[pos + len(ref):]        
        return sequence

    unalignedBam = tools.unalign_bam(bamFile, prefix + "_TEMP_unaligned")
    alignedBam = tools.align_pacbio(fa, unalignedBam, prefix + "_TEMP_pbsv", medianFilter=medianFilter)
    pbsvVCF = tools.find_sv(alignedBam, fa, prefix, indelOnly=True) 
    pbsvCallset = var.vcf_to_variantcallset(pbsvVCF)
    helper.delete_file(unalignedBam)
    
    seq = helper.get_fasta_seq(fa)
    
    seq = inject_sv(seq, pbsvCallset)
    fa = tools.dict2fasta(faDict, prefix + "TESTING")
    bam = tools.align_pacbio(fa, unalignedBam, prefix + "_TEMP_pbsv", medianFilter=False)
    pbsvVCF = tools.find_sv(bam, fa, prefix + "_TEST", indelOnly=True) 
    pbsvCallset = var.vcf_to_variantcallset(pbsvVCF)
    for variant in pbsvCallset.get_variant_calls(): print(variant)
    
    faDict = tools.fasta2dict(fa)
    fid = list(faDict.keys())[0]
    injectedFaDict = { fid + "|sv" : seq }
    injectedFa = tools.dict2fasta(injectedFaDict, prefix + "_pbsv_injected")
    
    return injectedFa



def inject_pbsv_calls_segment(segment, prefix, doA=True, doB=True):
    
    faA, faB = segment.write_seq(prefix)

    def inject_sv(sequence, callset):
        for pos in reversed(callset.get_positions()):
            variant = callset.get(pos)
            allele = variant.alleleB if variant.alleleA is None else variant.alleleA
            ref, alt = allele.get_ref(), allele.get_alt()
            sequence = sequence[:pos] + alt + sequence[pos + len(ref):]        
        return sequence

    if doA:
        unalignedBamA = tools.unalign_bam(segment.bamA, prefix + "_TEMP_unaligned_A")
        bamA = tools.align_pacbio(faA, unalignedBamA, prefix + "_TEMP_pbsv_A", medianFilter=False)
        pbsvVCF = tools.find_sv(bamA, faA, prefix + "_A", indelOnly=True) 
        pbsvCallset = var.vcf_to_variantcallset(pbsvVCF)
        helper.delete_file(bamA)
        
        seqA = segment.seqA
        
        seqA = inject_sv(seqA, pbsvCallset)
        faDict = {"a": seqA}
        fa = tools.dict2fasta(faDict, prefix + "TESTING")
        bam = tools.align_pacbio(fa, unalignedBamA, prefix + "_TEMP_pbsv_A", medianFilter=False)
        pbsvVCF = tools.find_sv(bam, fa, prefix + "_ATEST", indelOnly=True) 
        pbsvCallset = var.vcf_to_variantcallset(pbsvVCF)
        for variant in pbsvCallset.get_variant_calls(): print(variant)
        
    if doB:
        unalignedBamB = tools.unalign_bam(segment.bamB, prefix + "_TEMP_unaligned_B")
        bamB = tools.align_pacbio(faB, unalignedBamB, prefix + "_TEMP_pbsv_B", medianFilter=True)
        tools.find_sv(bamB, faB, prefix + "_B")
        pbsvVCF = tools.find_sv(bamB, faB, prefix + "_B", indelOnly=True)
        pbsvCallset = var.vcf_to_variantcallset(pbsvVCF)
        helper.delete_file(bamB)
        
        seqB = segment.seqB
        
        seqB = inject_sv(seqB, pbsvCallset)
        faDict = {"b": seqB}
        fa = tools.dict2fasta(faDict, prefix + "TESTING")
        bam = tools.align_pacbio(fa, unalignedBamB, prefix + "_TEMP_pbsv_B", medianFilter=True)
        pbsvVCF = tools.find_sv(bam, fa, prefix + "_BTEST", indelOnly=True) 
        pbsvCallset = var.vcf_to_variantcallset(pbsvVCF)
        print(pbsvCallset)

    
    return segment


def polish_junction(regionL, regionR, param, buffer=3000):
    
    bufferHalf = math.ceil(buffer/2)
    prefix = param.OUTPUT_DIR + "/polish_junction_TEMP_"   

    LA = regionL.seqA[-bufferHalf:] if regionL.seqA is not None else ""
    LAreads = tools.samtools_fetch(regionL.bamA) if regionL.bamA is not None else []
    LB = regionL.seqB[-bufferHalf:] if regionL.seqB is not None else ""
    LBreads = tools.samtools_fetch(regionL.bamB) if regionL.bamB is not None else []
    RA = regionR.seqA[:bufferHalf] if regionR.seqA is not None else ""
    RAreads = tools.samtools_fetch(regionR.bamA) if regionR.bamA is not None else []
    RB = regionR.seqB[:bufferHalf] if regionR.seqB is not None else ""
    RBreads = tools.samtools_fetch(regionR.bamB) if regionR.bamB is not None else []
    
    LUreads = tools.samtools_fetch(regionL.bamU) if regionL.bamU is not None else []
    RUreads = tools.samtools_fetch(regionR.bamU) if regionR.bamU is not None else []

    #todo: subset reads to the region we care about


    # left is unphased
    if len(LA) == 0 or len(LB) == 0:
        L = LB if len(LA) < 1 else LA
        A = L + RA
        B = L + RB
        
        Areads = LUreads + RAreads + RUreads
        Breads = LUreads + RBreads + RUreads
        
    # right is unphased
    elif len(RA) == 0 or len(RB) == 0:
        R = RB if len(RA) < 1 else RA
        A = LA + R
        B = LB + R
        
        Areads = LAreads + LUreads + RUreads
        Breads = LBreads + LUreads + RUreads
        
    # both are phased
    else:
        
        #todo: pair L and R in proper phase.
        A = LA + RA
        B = LB + RB 
        
        #todo: make read vector
        Areads = LAreads + RAreads # + UNPHASED
        Breads = LBreads + RBreads # + UNPHASED
    
    faDictA = {"A" : A}
    faDictB = {"B" : B}
    faA = tools.dict2fasta(faDictA, prefix + "A")
    faB = tools.dict2fasta(faDictB, prefix + "B")
    
    #todo: make sure headerbam exists
    unalignedBamA = tools.samtools_write(Areads, prefix + "A", regionL.bamU)
    bamA = tools.align_pacbio(faA, unalignedBamA, prefix + "A_aligned")    
    polishedFaA = iterative_polish(faA, bamA, prefix + "A", niter=2, bamUnphased=None)

    unalignedBamB = tools.samtools_write(Breads, prefix + "B", regionL.bamU)
    bamB = tools.align_pacbio(faB, unalignedBamB, prefix + "B_aligned")
    polishedFaB = iterative_polish(faB, bamB, prefix + "B", niter=2, bamUnphased=None)

'''

'''
    #todo: repetition causes expansion of junction because of the way reads align
    
    tools.fasta2dict(polishedFaB)["B|arrow|arrow"]
Out[1401]: 'CTATCCCAGACTTGCCCTCCTAAATGAAACCAAACTTATAAATGTTGTTAAAAGTTAAATTTATTTTAACTGATTGTTTTGCTTTCAAATACATATTATTGTTTTGCAAATTTATGTTATTAAAAATAAAATAAACAAGCCCAATGAAATTACCAAATTTTGATACATAGATTTGCCCAACTGGTTTAGCTTTAGATCCTCATCCAGTTGCTCAAGTGTTTTTTTGCAGTAAACAAGTTTCCTGTTTCAGCAGCTCCTGATGGTTGCTTTCAGCTTTCCTAGTTCTTTTTCATCTGTGATAGATTGAGCCGTATCTATGGCACATGACATAAGATTTCTCCTTTGACATCAAGTTCACTTCTAAGTTAATTCTGAAACTAACCTTCTGGAGGCAGTGTCAAAGTTTTTTTCTGAGAAGAGGGGTCAGTTCAATTGTACTATTAATTTTTAACTGTAATATATTTACTTGCATTTAGATTTTTCTTCAAGTTTGCATGATCTAGTTTTAATTGTCATTATGAGCAATTTCATCCCATTACTGTTTGTTTTTCTTGTACAAGTGATCATACCTTTAGTGTCATTTGAATTTTCAAAATAATAAAACCCATACTTTAAATAAGCATCAACATCCTACCTGTATTGACTACAATTGCTTGTAACCCATCTACCCATTATGGCTGAGGCAATATTCCCATTCCCCTTACTTTTGTTAATCCATATCTATCTTATGCCTTGCTGTTAGATTATTCTGAGGCGAAGAAAGAGAGGAAGAAGAAGAAAAGAAATAAGGAAGGAAGTGAGAAAGAGGGAAAGTGTTTCTCCATTCAAGGTATTATTATATCTAAAATTCACTAGTTGATTATCTAATATTTGTAAAATTACCTTCCAGATGAAGGAAATTATCTGATAAGATTCATGTAGGGGAAAGTCAAGGGAACTGTCATTAGGTGTCCATTAGATTCTTGGTACAATTTACCCTATTAATATGTGAGCTTGATTAAGGAGTTTTGGTGTCATATGTCAAAAAGTAACATATCTATTAAGAGAGAGGAAGGAAGTGGAGAGGGGAAAATATTGAGGGAAGGGTTGTTGAGATGGTTAAGAGAACTAGAGACTGATTCTTGATCTTTACTCAGGGCATATTTCTTGTTGAACAGCACTTTTCCCCTTCTGCTCATTCTCAAGGTGGTTTATCACTTCATAACTCCATGAATTCATTGCACAGGCCCATTGAATCACAAAACTATTTCTCACTTAAGTTCCCTACAGAAAAATCCCACAGCCACTGCAGTTTGGAAAGCAAAGGATTTAAATTATTTTTTCAAAAAGTGAGAATGAAAATAGTCCCATTGCACCCAGAAACCCTCCATGAATGTTGAGACAAGGCTAACTAAAAAGATTACCAGGGAGAGTAAAATGTGAGAGACGTGATTCCCCTTCCAGTAAATGCTTCCTTGGGAAAGGCGCTTGTTAACTCCCACTGAAAAATAGGCCCGCTGTCCAAGTCAGCAGTGAATTACAGAAGAAATGGTGCAAAACTTCACAGCCTGGAGAGAAAACCCCAGTTTCACTGAAGTGGAGCTGTGAAGTGATTCCCATATTTCAATAAGAATTATTTTATCTTTTTGCTAAGTCACTGAAATAAACAATCACATCTCTAAACTGTGAAAAACAGGAAGAAATACTGGGAAGATTTATAATACACATGTTTCTACTCCCTGTAGAAAAATGAAAATAGTAAATAAGTGCATTTCTGTAACATCCAAAATTAACAAAATAGTGTCTGTATGTTGTAGATAAAAAACCAATGTTCAAATAGGTTTAAGTGACTTCCCTATGGCCATAAAGGTGAGAAACGGTGATGTTGAAATGGCAACACTAGTCTTCTGAACTTCGGTCCATGCTTTTAAATCCCAGTCATAGCTGTAAGAGAAATATCTTAGAAGACAATTGTGTTATCATAATCCTGGATATGCTTATTTTTTTACCTGTATTATTTATTGGAGTTTTAAAAATGATTTTCCCATTTCTGTAATAAAAGTCAAAACATATTTGCTATAGTTTTATGTGTATTATTTTATAACATGCTATTCTCCCACTTAAGGATTTGGTTTGAGATCCAGATTTTGATATTCATTACCTCTAAGCAACCCTCTGTCATTTCTTACATAGAAGTTATGAATAATAACATTTGCCTTAGGATTATATTGTGTCTCTAGGTCCTAGTGAGTAGTGTGTGTTCAATAAATGGCAGATGTTATATTTGTTTCTATCCCAAAGCATGTGAATATACACTTTGTGCACATAAACAAAGGTTAATTAAAACAGAATTACATATATGCAAAATGAGAAACTAGAAGAAGGAGATGAGGAGGGTACTAAGAAACAACGGTAACAAGAAACACAATTAGCAAATCATTTTTTTTTTCAATCATGTGCAAGCCAGAAAAAAGGAAGTAACAATGTGGTTCTCATCAACTGGGTTAGGACCTCTTCTCTAGCAACCTGTATTTCTCTTTTTATAACACTCAGTTCATTTGATATATTTTTTCTTTTTGGTAATTTATTTAGGACTGTAAGCTCTATGAAGTCAAGCACTTTTTCCTTCCTGTATATCTAGCAGATAGAATAGAGTTGGGTACAGAGGACTCTCATCACATTTTGACAGCTCAAGTCTTTACTCCTCCCTCTCAAATAGCCCCCTGATTTGATTTTGGCAAATAATTCCTCTCTCCTTGTGTGATGTCTCAGTGGGAAAGTAAATCTATGTGCCAACATTTCACTATGTGGGTCCAAGTGATCAAATGCTTCTCTCTGGCTACAGTGGGGGCAAGCCTAATGCCTAAATTCAGCCAAACAAACACTCTCTCTTGAAAGTCTGAATTCTGAGTAAATGTCCACAACACAAAATAAATATTAGAATGTAATTATTTCACTTTCAGTGTCCTGATCAAACTATTGACTGCTATGCTATCATTGTTGTTGTTTCAAATAGGCATTG'

tools.fasta2dict(polishedFaB)["B|arrow|arrow"] == tools.fasta2dict(polishedFaA)["A|arrow|arrow"]
Out[1402]: True

L
Out[1403]: 'CTATCCCAGACTTGCCCTCCTAAATGAAACCAAACTTATAAATGTTGTTAAAAGTTAAATTTATTTTAACTGATTGTTTTGCTTTCAAATACATATTATTGTTTTGCAAATTTATGTTATTAAAAATAAAATAAACAAGCCCAATGAAATTACCAAATTTTGATACATAGATTTGCCCAACTGGTTTAGCTTTAGATCCTCATCCAGTTGCTCAAGTGTTTTTTTGCAGTAAACAAGTTTCCTGTTTCAGCAGCTCCTGATGGTTGCTTTCAGCTTTCCTAGTTCTTTTTCATCTGTGATAGATTGAGCCGTATCTATGGCACATGACATAAGATTTCTCCTTTGACATCAAGTTCACTTCTAAGTTAATTCTGAAACTAACCTTCTGGAGGCAGTGTCAAAGTTTTTTTCTGAGAAGAGGGGTCAGTTCAATTGTACTATTAATTTTTAACTGTAATATATTTACTTGCATTTAGATTTTTCTTCAAGTTTGCATGATCTAGTTTTAATTGTCATTATGAGCAATTTCATCCCATTACTGTTTGTTTTTCTTGTACAAGTGATCATACCTTTAGTGTCATTTGAATTTTCAAAATAATAAAACCCATACTTTAAATAAGCATCAACATCCTACCTGTATTGACTACAATTGCTTGTAACCCATCTACCCATTATGGCTGAGGCAATATTCCCATTCCCCTTACTTTTGTTAATCCATATCTATCTTATGCCTTGCTGTTAGATTATTCTGAGGCGAAGAAAGAGAGGAAGAAGAAGAAAAGAAATAAGGAAGGAAGTGAGAAAGAGGGAAAGTGTTTCTCCATTCAAGGTATTATTATATCTAAAATTCACTAGTTGATTATCTAATATTTGTAAAATTACCTTCCAGATGAAGGAAATTATCTGATAAGATTCATGTAGGGGAAAGTCAAGGGAACTGTCATTAGGTGTCCATTAGATTCTTGGTACAATTTACCCTATTAATATGTGAGCTTGATTAAGGAGTTTTGGTGTCATATGTCAAAAAGTAACATATCTATTAAGAGAGAGGAAGGAAGTGGAGAGGGGAAAATATTGAGGGAAGGGTTGTTGAGATGGTTAAGAGAACTAGAGACTGATTCTTGATCTTTACTCAGGGCATATTTCTTGTTGAACAGCACTTTTCCCCTTCTGCTCATTCTCAAGGTGGTTTATCACTTCATAACTCCATGAATTCATTGCACAGGCCCATTGAATCACAAAACTATTTCTCACTTAAGTTCCCTACAGAAAAATCCCACAGCCACTGCAGTTTGGAAAGCAAAGGATTTAAATTATTTTTTCAAAAAGTGAGAATGAAAATAGTCCCATTGCACCCAGAAACCCTCCATGAATGTTGAGACAAGGCTAACTAAAAAGATTACCAGGGAGAGTAAAATGTGAGAGACGTGATTCCCCTTCCAGTAAATGCTTCCTTGGGAAAGGCGCTTGTTAACTCCCACTGAAAAATAGGCCCGCTGT
    
    faDictLA = {str(regionL.regionA) + "A" : regionL.seqA}
    faDictLB = {str(regionL.regionB) + "B" : regionL.seqB}
    faDictRA = {str(regionR.regionA) + "A" : regionR.seqA}
    faDictRB = {str(regionR.regionB) + "B" : regionR.seqB}

    faLA = tools.dict2fasta(faDictLA, prefix + "LA")
    faLB = tools.dict2fasta(faDictLB, prefix + "LB")
    faRA = tools.dict2fasta(faDictRA, prefix + "RA")
    faRB = tools.dict2fasta(faDictRB, prefix + "RB")

    def align_fa(fa1, fa2):
        pafFile = tools.align_paf(fa1, fa2, prefix, secondary=False)
        return paf.get_variantset(pafFile, fa1)
        
    vsLALB = align_fa(faLA, faLB)
'''
    


