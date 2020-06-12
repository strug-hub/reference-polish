SAMTOOLS = "samtools"
PBMM2 = "pbmm2"
GENOMIC_CONSENSUS = "gcpp"
LONGSHOT = "longshot"
MM2 = "minimap2"
PAFTOOLS = "paftools"

import pysam
import pyfaidx
import copy
import subprocess
import os

def parse(command):
    return command.split()

#samtools
#========================================================

def samtools_index(bamFile):    
    """Index a bam file."""
    
    cmd = parse(SAMTOOLS)
    subprocess.call(cmd + ["index", bamFile])
    return bamFile + ".bai"

def samtools_fetch(bamFile, region=None, unaligned=True, asDict=False):
    samfile = pysam.AlignmentFile(bamFile, "rb")
    if region is None:
        alignments = [x for x in samfile.fetch(until_eof=unaligned)]
    else:
        alignments = [x for x in samfile.fetch(region.chrom, region.start, region.end)]
    
    if asDict:
        alignDict = dict()
        for a in alignments:
            if a.qname not in alignDict: alignDict[a.qname] = []
            alignDict[a.qname].append(a)
        return alignDict

    return alignments

def samtools_write(alignments, prefix, headerBam, makeUnique=False, index=True):
    outName = prefix + ".bam"

    # remove redundant reads:
    if makeUnique:
        readDict = {read.qname : read for read in alignments}
        uniqueReads = [readDict[k] for k in readDict]
        alignments = uniqueReads
        
    alignments = sorted(alignments, key=lambda x: (x.is_unmapped, x.rname, x.pos))

    header = pysam.AlignmentFile(headerBam, "rb")
    bamFile = pysam.AlignmentFile(outName, "wb", template=header)
    for alignment in alignments:
         bamFile.write(alignment)

    bamFile.close()
    header.close()
    
    if index:
        samtools_index(outName)
    return outName


#pacbio
#========================================================

def align_pacbio(refFa, readsFile, prefix, minConcordance=None, medianFilter=False):

    cmd = parse(PBMM2)
    outName = prefix + ".pbmm2.bam"
    
    pbmm2Align = cmd + ["align", "--sort"]
    
    if minConcordance is not None:
       pbmm2Align.extend(["--min-concordance-perc", str(minConcordance)])
       
    if medianFilter:
       pbmm2Align.append("--median-filter")

    pbmm2Align.append("--unmapped")

    pbmm2Align.extend([refFa, readsFile, outName])
    subprocess.call(pbmm2Align)
    samtools_index(outName)
    
    return outName

def pacbio_polish(bamFile, refFa, prefix, region=None, outputFasta=False, outputGff=False, useMapFilter=True, chunkSize=None):
    #samtools view -H $BAM | sed "s/VN:1.3/VN:1.3\tpb:3.0.4/" | samtools reheader - $BAM
        
    '''What is MapQV and why is it important?

    MapQV is a single scalar Phred-scaled QV per aligned read that reflects the mapper's degree of certainty that the read aligned to this part of the reference and not some other. Unambigously mapped reads will have a high MapQV (typically 255), while a read that was equally likely to have come from two parts of the reference would have a MapQV of 3.
    MapQV is pretty important when you want highly accurate variant calls. Quiver and Plurality both filter out aligned reads with a MapQV below 20 (by default), so as not to call a variant using data of uncertain genomic origin.
    This can be problematic if using quiver/arrow to get a consensus sequence. If the genome of interest contains long (relative to the library insert size) highly-similar repeats, the effective coverage (after MapQV filtering) may be reduced in the repeat regions---this is termed these MapQV dropouts. If the coverage is sufficiently reduced in these regions, quiver/arrow will not call consensus in these regions---see `What do quiver/arrow do for genomic regions with no effective coverage?`_.
    If you want to use ambiguously mapped reads in computing a consensus for a denovo assembly, the MapQV filter can be turned off entirely. In this case, the consensus for each instance of a genomic repeat will be calculated using reads that may actually be from other instances of the repeat, so the exact trustworthiness of the consensus in that region may be suspect. The next section describes how to disable the MapQV filter.
    How can the MapQV filter be turned off and when should it be?
    The MapQV filter can be disabled using the flag --mapQvThreshold=0 (shorthand: -m=0). If running a quiver/arrow job via SMRT Portal, this can be done by unchecking the "Use only unambiguously mapped reads" option. Consider this in de novo assembly projects, but it is not recommended for variant calling applications.
    '''    
    cmd = parse(GENOMIC_CONSENSUS)
    
    outVcf = prefix + ".consensus.vcf"
    outFa =  prefix + ".consensus.fasta"
    outGff = prefix + ".consensus.gff"

    outputFiles = outVcf

    if outputFasta:
        outputFiles = outputFiles + "," + outFa
    if outputGff:
        outputFiles = outputFiles + "," + outGff
    arrow = cmd + ["-r", refFa, "-o", outputFiles ]
    arrow.extend(["--algorithm", "arrow"])
    if region is not None:
        arrow.extend(["--windows", str(region)])
        
    if not useMapFilter:
        arrow.append("-m=0")

    if chunkSize is not None:
        arrow.extend(["-C", str(chunkSize)])

    arrow.append(bamFile)
    
    print(arrow)
    subprocess.call(arrow)
    
    if outputFasta: 
        pyfaidx.Faidx(outFa)    
        return outFa
    
    return outVcf

#longshot
#========================================================

LONGSHOT_BAM_SUFFIX = ".longshot.bam"

def longshot_genotype(bamFile, refFa, prefix, region=None, coverageAware=True, writeBams=False):
    
    outName = prefix + ".longshot.vcf"
    cmd = parse(LONGSHOT)
    longshot = cmd
    
    if coverageAware:
        longshot.append("-A")
    
    if writeBams:
        #old version of longshot
        #longshot.extend(["--hap_bam_prefix", prefix + LONGSHOT_BAM_PREFIX])
        longshot.extend(["--out_bam", prefix + LONGSHOT_BAM_SUFFIX])

    #force overwrite of output file
    longshot.append("-F")

    if region is not None:
        longshot.extend(["--region", region.str_base1()])

    longshot.extend(["--bam", bamFile, "--ref", refFa, "--out", outName])
    subprocess.call(longshot)
    
    
    if writeBams:
        '''
        bamFilePrefix = prefix + LONGSHOT_BAM_PREFIX
        hpA = bamFilePrefix + ".hap1.bam"
        hpB = bamFilePrefix + ".hap2.bam"
        unphased = bamFilePrefix + ".unassigned.bam"
        samtools_index(hpA)
        samtools_index(hpB)
        samtools_index(unphased)
        '''
        samtools_index(prefix + LONGSHOT_BAM_SUFFIX)

    return outName

def get_longshot_phased_reads(prefix, keepOriginal=False):
   
    #old version of longshot
    #--hap_bam_prefix
    #Write haplotype-separated reads to 3 bam files using this prefix:
    #<prefix>.hap1.bam, <prefix>.hap2.bam, <prefix>.unassigned.bam

    '''
    bamFilePrefix = prefix + LONGSHOT_BAM_PREFIX
    h1= bamFilePrefix + ".hap1.bam"
    h2 = bamFilePrefix + ".hap2.bam"
    unphased = bamFilePrefix + ".unassigned.bam"
    '''

    h1 = prefix + ".longshot.hap1"
    h2 = prefix + ".longshot.hap2"
    unphased = prefix + ".longshot.unassigned"

    if os.path.isfile(h1 +".bam" ) and os.path.isfile(h2 +".bam" ) and \
        os.path.isfile(unphased +".bam" ):
            return [x+".bam" for x in (h1, h2, unphased)]

    bamFile = prefix + LONGSHOT_BAM_SUFFIX 
    hap = {x:[] for x in ["0","1","2"]}

    
    for read in samtools_fetch(bamFile):
        try:
            hp = str(read.get_tag("HP"))
        except:
            hp = "0"
        
        hap[hp].append(read)
        

    samtools_write(hap["1"], h1, bamFile)
    samtools_write(hap["2"], h2, bamFile)
    samtools_write(hap["0"], unphased, bamFile)
    
    if not keepOriginal:
        try:
            os.remove(bamFile)
            os.remove(bamFile + ".bai")
        except:
            pass
        
    return [x+".bam" for x in (h1, h2, unphased)]

#minimap2
#========================================================

def align_paf(refFa, readsFile, prefix, asm="asm5", r=500, secondary=False):
    '''
    -r
    Bandwidth used in chaining and DP-based alignment [500].
    This option approximately controls the maximum gap size.
    '''
    cmd = parse(MM2)
    outName = prefix + ".paf"
    
    mm2Align = cmd + ["-cx", asm, "--cs", "-r", str(r)]
    
    if not secondary:
        mm2Align.append("--secondary=no")
    
    mm2Align.extend(["-o", outName, refFa, readsFile])
    subprocess.call(mm2Align)
    
    return outName

def align_paf_lenient(refFa, readsFile, prefix):
    return align_paf(refFa, readsFile, prefix, r=10000, asm="asm10", secondary=False)
def align_paf_very_lenient(refFa, readsFile, prefix):
    return align_paf(refFa, readsFile, prefix, r=20000, asm="asm20", secondary=False)

def paf_liftover(pafFile, regions, minMapQuality=5, maxDivergance=1, minAlnLen=2500):
    '''
    -q INT    min mapping quality [5]
    -l INT    min alignment length [50000]
    -d FLOAT  max sequence divergence (>=1 to disable) [1]
    '''
    cmd = parse(PAFTOOLS)
    bedFile = os.path.dirname(os.path.realpath(pafFile)) + "/_tempbed_.bed"
    
    if not isinstance(regions, list): regions = [regions]

    f = open(bedFile, "w")
    for region in regions:
        f.write("\t".join([str(x) for x in (region.chrom, region.start, region.end)]) + "\n")
    f.close()
    
    paftoolsLiftover = cmd + ["liftover"]  + \
                    ["-d", str(maxDivergance)] + \
                    ["-l", str(minAlnLen)] + \
                    ["-q", str(minMapQuality)]

    paftoolsLiftover.extend([pafFile])
    paftoolsLiftover.extend([bedFile])
    
    outFile = os.path.dirname(os.path.realpath(pafFile)) + "/_result_.bed"

    writer = open(outFile, 'w+')
    subprocess.call(paftoolsLiftover, stdout=writer)
    writer.close()

    liftover = []    
    for line in open(outFile, "r"):
        r = copy.deepcopy(regions[0])
        s = line.split("\t")
        r.chrom = s[0] ; r.start = int(s[1]); r.end = int(s[2])
        liftover.append(r)


    #todo: outfile missing matches?
    os.remove(bedFile)
    os.remove(outFile)

    return liftover

def paf_call(pafFile, referenceFa, prefix, vcfSample="sample", minAlnLen=2500):
    '''
      -f FILE   reference sequences (enabling VCF output) [null]
      -s NAME   sample name in VCF header [sample]
    '''
    cmd = parse(PAFTOOLS)

    paftoolsCall = cmd + ["call"]  + \
                    ["-f", str(referenceFa)] + \
                    ["-L", str(minAlnLen)] + \
                    ["-s", str(vcfSample)]

    paftoolsCall.extend([pafFile])
    
    outFile = prefix + ".vcf"

    print(" ".join(paftoolsCall) )
    writer = open(outFile, 'w+')
    subprocess.call(paftoolsCall, stdout=writer)
    writer.close()

    return outFile
