import re
import random, math

import utils.external_tools as tools
import utils.file_handler as io
import utils.fasta_handler as fasta
from utils.region import SimpleRegion

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


