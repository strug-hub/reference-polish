import vcf
import os
import re
import random
import math
import sys
sys.path.append("..")

import utils.file_handler as io
import utils.fasta_handler as fasta
import utils.vcf_annotation as annotator
from structures.region import SimpleRegion
import utils.external_tools as tools
import analysis.analyze_sequence as analyzer

def fetch_ref_reads(bamFile, region, outdir, fileSuffix=""):
    refAlignments = tools.samtools_fetch(bamFile, region)
    return tools.samtools_write(refAlignments, outdir + "pacbio_reads" + fileSuffix, bamFile)

def fetch_query_reads(bamFile, region, outdir):
    io.delete_file(outdir + "10x_reads", deleteDirTree=True)
    return tools.bam2fastq_10x(bamFile, outdir, dirName="10x_reads", region=region, returnFqList=True)

def align_ref(faFile, refReads, outdir, param, outName="ref_aligned"):
    
    bam = outdir + outName + ".pbmm2.bam"
    if os.path.isfile(bam):
        print("pbmm2 alignments found, skipping step")
        return bam
    
    bam = tools.align_pacbio(faFile, refReads, outdir + outName)
    return bam

def align_query(faFile, queryReads, outdir, param, outName="query_aligned"):
    #note: outName is irrelevant
    
    queryReadsDir = "/".join(queryReads[0].split("/")[:-1])
    
    sampleName = ".".join(faFile.split(".")[:-1]).split("/")[-1] + "-longranger"

    longrangerBam = outdir + sampleName + "/phased_possorted_bam.bam"
    longrangerVCF = outdir + sampleName + "/phased_variants.vcf.gz"

    if not os.path.isfile(longrangerBam) or not os.path.isfile(longrangerVCF):
        longrangerVCF, longrangerBam = tools.align_10x(faFile, queryReadsDir, outdir, svBlacklist=True)
    else:
        print("Long Ranger files found, skipping step")
        return longrangerBam

    #correct bug in longranger 2.2.2 where hets are sometimes given GT = 1|1 

    correctedVCF = re.sub(".gz", "", longrangerVCF)
    annotator.run_correction(longrangerVCF, correctedVCF)
    io.delete_file(longrangerVCF)
    
    longrangerVCF = tools.bgzip(correctedVCF)
    #longrangerVariants = [record for record in vcf.Reader(open(longrangerVCF, "rb"))]

    if not param.KEEP_INTERMEDIATE:
        io.delete_file(outdir + "refdata-" + outName, deleteDirTree=True)

    return longrangerBam


def polish_ref(targetFa, region, outdir, param, refReads=None, refAlignments=None, outName="ref_polished"):
        
    polishedFa = outdir + outName + ".consensus.fasta"
    if os.path.isfile(polishedFa):
        print("Arrow polished fasta found, skipping step")
        return polishedFa
    
    deleteAlignments = False
    if refAlignments is None:
        refAlignments = align_ref(targetFa, refReads, outdir, param, outName)
        deleteAlignments = True

    polishedFa = tools.pacbio_polish(refAlignments, targetFa, outdir + outName, region=region, outputFasta=True)

    if deleteAlignments:
        io.delete_file(refAlignments)
        
    fasta.index_fasta(polishedFa)
    return polishedFa

def polish_query(targetFa, region, outdir, param, queryReads=None, queryAlignments=None, outName="query_polished"):

    polishedFa = outdir + outName + ".fasta"
    if os.path.isfile(polishedFa):
        print("Pilon polished fasta found, skipping step")
        return polishedFa

    deleteAlignments = False
    if queryAlignments is None:
        qreads1, qreads2 = queryReads
        trimmedr1 = tools.trim_10x_barcode(qreads1)
        queryAlignments = tools.bwa_mem(targetFa, trimmedr1, qreads2, outdir + "bwa_aligned")
        deleteAlignments = True

    pilonFa = tools.pilon_polish(queryAlignments, targetFa, outdir, outputVCF=False)
    polishedFa = outdir + outName + ".fasta"

    io.move_file(pilonFa, polishedFa)
    fasta.index_fasta(polishedFa)
    
    for ext in [".sa", ".pac", ".bwt", ".ann", ".amb", ".ann"]:
        io.delete_file(polishedFa + ext)
        
    if deleteAlignments:
        io.delete_file(queryAlignments)      
        io.delete_file(trimmedr1)
        
    return polishedFa

def high_conf_hets(consensusFa, refBam, queryBam, outdir, param):
    
    vcfFile = outdir + "high_confidence_hets.vcf"
    if os.path.isfile(vcfFile):
        print("High confidence heterozygous variants found, skipping step")
        return vcfFile

    longshotVCF = phase_region(consensusFa, refBam, outdir + "phased", writeBams=False)
    longrangerVCF = "/".join(queryBam.split("/")[:-1]) +  "/phased_variants.vcf.gz" 

    longshotVariants = [record for record in vcf.Reader(open(longshotVCF, "rb" if longshotVCF.endswith("gz") else "r"))]
    longrangerVariants = [record for record in vcf.Reader(open(longrangerVCF, "rb" if longrangerVCF.endswith("gz") else "r"))]

    longshotSnps = {v.POS : v for v in longshotVariants if v.is_snp and v.samples[0].is_het}
    longrangerSnps = {v.POS : v for v in longrangerVariants if v.is_snp and v.samples[0].is_het}
    
    commonSnpPos = set()
    commonSnps = []

    for pos in longshotSnps:
        if pos in longrangerSnps and \
           longshotSnps[pos].REF == longrangerSnps[pos].REF and \
           longshotSnps[pos].ALT == longrangerSnps[pos].ALT:
               
               if len(longshotSnps[pos].FILTER) < 1 and len(longrangerSnps[pos].FILTER) < 1:
                   commonSnpPos.add(pos)
                   commonSnps.append(longrangerSnps[pos])

    longrangerOnly = [longrangerSnps[pos] for pos in longrangerSnps if pos not in commonSnpPos]
    highQualityLongranger = [v for v in longrangerOnly if len(v.FILTER) < 1 and v.QUAL > 150]

    longshotOnly = [longshotSnps[pos] for pos in longshotSnps if pos not in commonSnpPos]
    
    #TODO: find a good threshold here for quality
    highQualityLongshot = [v for v in longshotOnly if len(v.FILTER) < 1 and v.QUAL > 15]

    highConfidenceSnps = commonSnps + highQualityLongranger + highQualityLongshot
    
    longrangerNotSnp = [v for v in longrangerVariants if not v.is_snp and v.samples[0].is_het]
    
    highConfidenceIndels=[]
    for v in longrangerNotSnp:
        
        #bialleleic
        if len(v.alleles) > 2 or len(v.ALT) > 1: 
            continue
        #homopolymer
        if v.INFO["POSTHPC"][0] > 1:
            continue
        #filtered
        if len(v.FILTER) > 0:
            continue
        #complex
        if len(v.REF) != 1 and len(v.ALT[0]) != 1:
            continue
        if len(v.REF) > 8 or len(v.ALT[0]) > 8:
            continue
        #low qual
        if v.QUAL < 100:
            continue

        if not v.samples[0].is_het:
            continue
                
        highConfidenceIndels.append(v)

    highConfidenceVariants = highConfidenceSnps + highConfidenceIndels
    
    #vcfFile = outdir + "high_confidence_hets.vcf"
    writer = open(vcfFile, "w+")
    writer.write("##fileformat=VCFv4.2\n")
    header = ["#CHROM", "POS", "ID", "REF", "ALT",
               "QUAL",  "FILTER", "INFO", "FORMAT", "SAMPLE"]
    writer.write("\t".join(header) + "\n")

    for variant in sorted(highConfidenceVariants, key=lambda x: x.POS):
        line = [variant.CHROM, variant.POS, ".", variant.REF, variant.ALT[0], 
                variant.QUAL, "PASS", ".", "GT", "0/1"]       
        writer.write("\t".join([str(x) for x in line]) + "\n")

    writer.close()
    
    return vcfFile

def phase_consensus(consensusFa, highConfVCF, refBam, queryBam, outdir, param):

    #phase het variants with all reads using whatshap
    whatshapVCF = tools.whatshap_phase([refBam, queryBam], highConfVCF, consensusFa,
                                 outdir + "high_confidence",  indels=True, maxCov=17)
    
    refHaploBam = tools.whatshap_haplotag(refBam, whatshapVCF, consensusFa)    
    queryHaploBam = tools.whatshap_haplotag(queryBam, whatshapVCF, consensusFa)
    
    if not param.KEEP_INTERMEDIATE:
        io.delete_file(refBam)
    
    return refHaploBam, queryHaploBam        
    

def haplotype_polish_query(hap, hapFa, haploBam, outdir, param, realign=True):
    hap = str(hap)
    reads = split_reads(haploBam)

    assignedReads = reads[hap]
    idSet = {r.query_name for r in assignedReads}

    for r in reads["1" if hap == "2" else "2"] + reads["0"]:
        if r.query_name in idSet:
            assignedReads.append(r)
    
    assignedBam = tools.samtools_write(assignedReads, outdir + hap + "_polish_input_query", haploBam)
    
    if realign:
        dirname = "10x_hap" + hap
        r1,r2 = tools.bam2fastq_10x(assignedBam, outdir, dirName=dirname, returnFqList=True)
        r1Trimmed = tools.trim_10x_barcode(r1)
        tools.bwa_index(hapFa)
        realignedBam = tools.bwa_mem(hapFa, r1Trimmed, r2, outdir + "bwa_aligned_hap" + hap)
        
        io.delete_file(outdir+dirname, deleteDirTree=True)
        io.delete_file(assignedBam)
        assignedBam = realignedBam
        for ext in [".sa", ".pac", ".bwt", ".ann", ".amb", ".ann"]:
            io.delete_file(hapFa + ext)

    hapPolishedFa = tools.pilon_polish(assignedBam, hapFa, outdir, prefix="hap" + hap + ".consensus", 
                                       changeFile=True, snpIndelOnly=True)
    fasta.index_fasta(hapPolishedFa)

    #helper.delete_file(assignedBam)
    return hapPolishedFa

def haplotype_polish_ref(hap, hapFa, haploBam, outdir, param, realign=True):

    hap = str(hap)
    reads = split_reads(haploBam)

    random.shuffle(reads["0"])
    half = math.ceil(len(reads["0"])/2)
    assignedReads = reads[hap] + reads["0"][:half]
    assignedBam = tools.samtools_write(assignedReads, outdir + hap + "_polish_input_ref", haploBam)

    if realign:
        realignedBam = tools.align_pacbio(hapFa, assignedBam, outdir + hap + "_polish_input_ref.realigned")
        io.delete_file(assignedBam)
        assignedBam = realignedBam
        
    hapPolishedFa = tools.pacbio_polish(assignedBam, hapFa, outdir + "hap" + hap, outputFasta=True)
    
    io.delete_file(assignedBam)

    return hapPolishedFa

def clean_directory(outdir):
    
    #delete reads
    io.delete_file(outdir + "10x_reads", deleteDirTree=True)
    io.delete_file(outdir + "pacbio_reads.bam")
    

def reference_analysis(faFile, refReads, outdir, param):
    
    outName="analysis"
    refBam = tools.align_pacbio(faFile, refReads, outdir + outName)
    alignments = tools.samtools_fetch(refBam)


    alignments = tools.samtools_fetch(param.REF_ALIGNED_READS, rRegion)
    analyzer

    return refBam, queryBam

def analyze_contig(refReads, queryReads, unpolishedFa, consensusFa, hap1Fa, hap2Fa, outDir, param):
    
    d="/media/scott/Zapdos/out/scaff139_/"
    outDir=d
    consensusFa=outDir+"consensus.fasta"     
    unpolishedFa=outDir+"unpolished.fasta" 
    consensusFa=outDir+"consensus.fasta" 
    hap1Fa=outDir+"hap1.fasta" 
    hap2Fa=outDir+"hap2.fasta" 

    refReads=outDir+"pacbio_reads.bam" 
    queryReads=(outDir+"10x_reads/bamtofastq_S1_L000_R1_001.fastq.gz", outDir+"10x_reads/bamtofastq_S1_L000_R2_001.fastq.gz")
    faFile=unpolishedFa

    phasedBam=outDir+"ref_to_consensus.pbmm2.haplotag.bam"
    phasedAlns = tools.samtools_fetch(phasedBam)
    h1,h2,unknown = [],[],[]
    for aln in phasedAlns:
        try:
            hp=str(aln.get_tag("HP"))
            if hp == "1":
                h1.append(aln)
            elif hp == "2":
                h2.append(aln)
            else:
                unknown.append(aln)
        except:
            unknown.append(aln)

    tempDir = outDir + ("/" if outDir[-1] != "/" else "") + "temp/"
    io.make_dir(tempDir)

    h1RefBam = tools.samtools_write(h1, tempDir + "h1_ref", phasedBam, index=True)
    h2RefBam = tools.samtools_write(h2, tempDir + "h2_ref", phasedBam, index=True)
    h0RefBam = tools.samtools_write(unknown, tempDir + "h0_ref", phasedBam, index=True)

    files = [unpolishedFa, consensusFa, hap1Fa, hap2Fa]
    names = ["Unpolished", "Consensus", "Hap1", "Hap2"]
    
    refBamFiles, queryBamFiles = [], []
    
    concordanceList, concordanceListName = [], []
    for faFile, name in zip(files, names):
        haplotypes = ["1","2","0"]
        bamFiles = []
        for haplotypeBam, hp in zip([h1RefBam, h2RefBam, h0RefBam], haplotypes):
    
            if not os.path.isfile(tempDir + name + "." + hp + ".pbmm2.bam"):
                rBam = tools.align_pacbio(faFile, haplotypeBam, tempDir + name + "." +hp)
                #analyzer.plot_coverage(rBam, title=None)
            else:
                rBam = tempDir + name + "." + hp + ".pbmm2.bam"
            refBamFiles.append(rBam)
            concordanceList.append(analyzer.get_concordance(rBam, hap=hp))
            concordanceListName.append(name)
            bamFiles.append(rBam)
            
        #analyzer.get_coverage(bamFiles, haplotypes, outDir)

            
            
    analyzer.plot_concordance(concordanceList, concordanceListName, outDir)




    
    '''
        qreads1, qreads2 = queryReads
        trimmedr1 = tools.trim_10x_barcode(qreads1)
        queryBam = tools.bwa_mem(faFile, trimmedr1, qreads2, tempDir + "bwa_aligned" + name)
        io.delete_file(trimmedr1)
        
        queryBamFiles.append(queryBam)
    '''

        

            
            
            
            
            
    

        

    for faFile, name in zip(files, names):
        refBam, queryBam = align_for_analysis(faFile, h1Bam, queryReads, outDir, param, outName=name)
        concordanceList.append(analyzer.get_concordance(refBam))
        #analyzer.plot_coverage(queryBam, title=None)



    analyzer.plot_concordance(concordanceList, names, title=None)


def polish_across_gap(gapFa, refBam, seqData, outdir, param, outName="ref_polished"):
        
    polishedFa = tools.pacbio_polish(refBam, gapFa, outdir + outName + "_temp", outputFasta=True)
    io.delete_file(re.sub(".fasta", ".vcf", polishedFa))       
    
    return polishedFa

def align_across_gap_ref(gapFa, refReads, rRegion, seqData, outdir, param, fileSuffix=""):
    gapBam = tools.align_pacbio(gapFa, refReads, outdir + "pacbio_aligned_gap_" + fileSuffix)
    return gapBam

def phase_vcf(consensusFa, consensusVariants, refBam, queryBam, outdir, param):
    
    filteredVCF = annotator.filter_nonpass(consensusVariants, outdir + "variant_calls.filtered.vcf")
    whatshapVCF = tools.whatshap_phase([refBam, queryBam], filteredVCF, consensusFa,
                                 outdir + "variant_calls.filtered.phased", genotype=False, indels=True, maxCov=17)
    
    return whatshapVCF


def cluster_variants(variants1, variants2, window):

    vcollapse1 = []
    startNewCollapseGroup = True
    group = []
    for variant in variants1:
        
        if len(group) > 0 and abs(variant.POS - group[-1].POS) <= window:
            group.append(variant)
        else:
            startNewCollapseGroup = True

        if startNewCollapseGroup:
            if len(group) > 0: vcollapse1.append(group)
            group = [variant]
            startNewCollapseGroup = False

    if len(group) > 0: vcollapse1.append(group)

    vcollapse2 = []
    for group in vcollapse1:
        start, end = group[0].POS, group[-1].POS
        
        paired = []
        
        i =0
        while i < len(variants2):
            if variants2[i].POS >= (start-window) and variants2[i].POS <= (end+window):
                paired.append(variants2.pop(i))
                continue
            i=i+1
        
        vcollapse2.append(paired)

    startNewCollapseGroup = True
    group = []
    for variant in variants2:
        
        if len(group) > 0 and abs(variant.POS - group[-1].POS) <= window:
            group.append(variant)
        else:
            startNewCollapseGroup = True

        if startNewCollapseGroup:
            if len(group) > 0: 
                vcollapse1.append([])              
                vcollapse2.append(group)

            group = [variant]
            startNewCollapseGroup = False
    if len(group) > 0: 
        vcollapse2.append(group)
        vcollapse1.append([])              

    return (vcollapse1, vcollapse2)


def rescue_variant(vcf1, vcf2, mappableRegions):
    
    variants1 = [record for record in vcf.Reader(open(vcf1, "rb" if vcf1.endswith("gz") else "r"))]
    variants1 = sorted(variants1, key=lambda x: x.POS)
    variants2 = [record for record in vcf.Reader(open(vcf2, "rb" if vcf2.endswith("gz") else "r"))]
    variants2 = sorted(variants2, key=lambda x: x.POS)

    # replace the None caused by "." with an empty list
    for v in variants1 + variants2:
        if v.FILTER is None: v.FILTER = []
    
    def in_mappable_region(v):
        for mr in mappableRegions:
            if mr.contains_pos(v.POS):
                return True
        return False

    isMappable1 = [in_mappable_region(v) for v in variants1]
    isMappable2 = [in_mappable_region(v) for v in variants2]
    
    #rescue pacbio variants in difficult regions
    rescue1 = [v for mappable,v in zip(isMappable1, variants1) if not mappable]
    variants1 = [v for mappable,v in zip(isMappable1, variants1) if mappable]
    #filter out pacbio homopolymers
    filter1 = [v for v in variants1 if "HOMP" in v.INFO]
    variants1 = [v for v in variants1 if "HOMP" not in v.INFO]

    #filter 10x variants in difficult regions
    filter2 = [v for mappable,v in zip(isMappable2, variants2) if not mappable]
    variants2 = [v for mappable,v in zip(isMappable2, variants2) if mappable]
    #filter out 10x flagged variants
    filter2 = [v for v in variants2 if v.FILTER is not None]
    variants2 = [v for v in variants2 if len(v.FILTER) == 0 ]
    
    #rescue 10x homopolymers unless they are very long
    MAX_HOMOPOLYMER_LENGTH=50
    rescue2, leftover2 = [],[]
    for v in variants2:
        if "HOMP" in v.INFO and v.INFO["HOMP"] + v.INFO["HOMP_DIFF"] <= MAX_HOMOPOLYMER_LENGTH:
            rescue2.append(v)
        else: leftover2.append(v)

    variants2 = leftover2


    vcollapse1, vcollapse2 = cluster_variants(variants1, variants2, window=121)

    def get_change(variants):
        insDict = {"A":0, "T":0, "G":0, "C":0}
        delDict = {"A":0, "T":0, "G":0, "C":0}

        for v in variants:
            ref, alt=str(v.REF), str(v.ALT[0])
            while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
                ref, alt = ref[1:], alt[1:]            
            while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
                ref, alt = ref[:-1], alt[:-1]
            for char in ref:
                if char in delDict: delDict[char] +=1
            for char in alt:
                if char in insDict: insDict[char] +=1

        return(insDict, delDict)
    
    variants1, variants2 = [],[]
    for g1, g2 in zip(vcollapse1, vcollapse2):
        
        if len(g1) > 0 and len(g2) > 0:

            '''
            print("-----------------------------")
            print("\n".join([str(g) + "\t" + str(g.samples[0].gt_alleles) for g in g1]), "\n")
        
            print("***")
            print("\n".join([str(g) + "\t" + str(g.samples[0].gt_alleles) for g in g2]), "\n")
            '''
            
            ins1, del1 = get_change(g1)
            change1 = sum([ins1[x] + del1[x] for x in ["A","T","C","G"]])
            ins2, del2 = get_change(g2)
            change2 = sum([ins2[x] + del2[x] for x in ["A","T","C","G"]])
            
            editDist = sum([abs(ins2[x] - ins1[x]) + abs(del2[x] - del1[x]) for x in ["A","T","C","G"]])
            netEffect = editDist - max(change1, change2)
            
            #todo: recluster w smaller window?

            #some variation detected in both, trust 10x in this cluster
            if netEffect < 1:
                rescue2.extend(g2)
                filter1.extend(g1)
                continue
            
        filter1.extend(g1)
        filter2.extend(g2)

    return (rescue1, filter1, rescue2, filter2)
    

def call_variants(consensusFa, queryCallsVCF, refCallsVCF, mappableRegions, outdir):
        
    #happyVCF = tools.hap_py(refFa, longrangerVCF, mergedVCF, prefix, engine="vcfeval")
    rtgVCF = tools.vcfeval(consensusFa, queryCallsVCF, refCallsVCF, outdir)
    #tp = rtgVCF + "tp.vcf.gz"
    tpBL = rtgVCF + "tp-baseline.vcf.gz"
    fn = rtgVCF + "fn.vcf.gz"
    fp = rtgVCF + "fp.vcf.gz"
    
    intersectingVariants = [record for record in vcf.Reader(open(tpBL, "rb" if tpBL.endswith("gz") else "r"))]

    (rescue1, filter1, rescue2, filter2) = rescue_variant(fp, fn, mappableRegions)

    passVariants = intersectingVariants + rescue1 + rescue2
    failedVariants = filter1 + filter2
    
    for variant in passVariants:
        variant.FILTER = []
    for variant in failedVariants:
        variant.FILTER = ["FAILED"]

    outVCF = outdir + "variant_calls.vcf"
    header = vcf.Reader(open(queryCallsVCF, "rb" if queryCallsVCF.endswith("gz") else "r"))
    writer = vcf.Writer(open(outVCF, "w"), header)

    for variant in passVariants+failedVariants:
        writer.write_record(variant)
    
    writer.flush()
    writer.close()

    return outVCF

def call_variants_query(consensusFa, queryBam, outdir, param):
    longrangerVCF= "/".join(queryBam.split("/")[:-1]) + "/phased_variants.vcf.gz" 
    
    normalizedLongranger = tools.pre_py(consensusFa, longrangerVCF, noGz=True)
    annotatedLongrangerVCF = annotator.annotate_hompolymers(normalizedLongranger, outdir + "query_phased.vcf", consensusFa)
    annotatedLongrangerVCF = tools.bgzip(annotatedLongrangerVCF)
    
    io.delete_file(normalizedLongranger)
    return annotatedLongrangerVCF

def split_reads(bamFile):
    reads = tools.samtools_fetch(bamFile)
    readDict = {"1" : [], "2" : [], "0" : []}
    
    for read in reads:
        try:
            hp = read.get_tag("HP")
            if hp == 1: readDict["1"].append(read)
            elif hp == 2: readDict["2"].append(read)
            else: readDict["0"].append(read)
        except:
            readDict["0"].append(read)
    return readDict


def call_variants_ref(consensusFa, refBam, outdir, param):
        
    reads = split_reads(refBam)

    def haplo_polish(hap, reads):
        hapPolishedFa = outdir + hap + "_polish.consensus.fasta"
        hapPolishedVCF = outdir + hap + "_polish.consensus.vcf"
        if os.path.isfile(hapPolishedFa) and os.path.isfile(hapPolishedVCF):
            print("Haplotype " + hap + " polished files found, skipping step")
            return [hapPolishedFa, hapPolishedVCF]
        
        random.shuffle(reads["U"])
        half = math.ceil(len(reads["U"])/2)
        assignedReads = reads[hap] + reads["U"][:half]
        assignedBam = tools.samtools_write(assignedReads, outdir + hap + "_polishinput", refBam)
        hapPolishedFa = tools.pacbio_polish(assignedBam, consensusFa, outdir + hap + "_polish", outputFasta=True)
        io.delete_file(assignedBam)
        hapPolishedVCF = re.sub(".fasta", ".vcf", hapPolishedFa)
        return [hapPolishedFa, hapPolishedVCF]
    

    fastaA,vcfA = haplo_polish("A", reads)
    fastaB,vcfB = haplo_polish("B", reads)

    #optional back align step
    tools.align_pacbio(consensusFa, fastaA, outdir + "hapA_backaligned")
    tools.align_pacbio(consensusFa, fastaB, outdir + "hapB_backaligned")

    #annotate vcfs
    vcfA_GT = annotator.add_gt(vcfA, outdir + "A.gt.vcf")
    vcfB_GT = annotator.add_gt(vcfB, outdir + "B.gt.vcf")
    normalizedA = tools.pre_py(consensusFa, vcfA_GT, noGz=True)
    normalizedB = tools.pre_py(consensusFa, vcfB_GT, noGz=True)
    
    tempVCF = annotator.combine_vcf_lines(normalizedA, normalizedB, outdir + "merged_temp.vcf")
    mergedVCF = annotator.annotate_hompolymers(tempVCF, outdir + "ref_phased.vcf", consensusFa)
    mergedVCF = tools.bgzip(mergedVCF)

    io.delete_file(vcfA_GT) ;     io.delete_file(vcfB_GT)
    io.delete_file(normalizedA) ; io.delete_file(normalizedB)
    io.delete_file(tempVCF)

    return mergedVCF
   



def phase_reads(consensusFa, highConfVCF, refBam, queryBam, outdir, param):
    
    whatshapVCF = tools.whatshap_phase([refBam, queryBam], highConfVCF, consensusFa,
                                 outdir + "high_confidence",  indels=True, maxCov=17)

    haplotaggedBam = tools.whatshap_haplotag(refBam, whatshapVCF, consensusFa)
    
    return haplotaggedBam, queryBam


def phase_region(fa, reads, prefix, region=None, writeBams=True):
    longshotVCF = tools.longshot_genotype(reads, fa, prefix, region, writeBams=writeBams)
    if not writeBams: return longshotVCF
    bamA, bamB, bamUnphased = tools.get_longshot_phased_reads(prefix)
    return (longshotVCF, bamA, bamB, bamUnphased)

def phase_reads_refonly(consensusFa, refBam, outdir, param):
    
    longshotVCF, bamA, bamB, bamUnphased = phase_region(consensusFa, refBam, outdir + "phased", writeBams=True)
    
    A = tools.samtools_fetch(bamA)
    B = tools.samtools_fetch(bamB)
    U = tools.samtools_fetch(bamUnphased)

    for read in A: read.set_tag("HP", 1, "i")
    for read in B: read.set_tag("HP", 2, "i")
    
    haplotaggedBam = tools.samtools_write(A+B+U, re.sub(".bam", ".haplotag", refBam), refBam)
    return haplotaggedBam




def align_to_consensus_query(consensusFa, queryReads, qRegion, seqData, outdir, param):
    queryReadsDir = "/".join(queryReads[0].split("/")[:-1])
    
    longrangerBam = outdir + "consensus-longranger/phased_possorted_bam.bam"
    longrangerVCF = outdir + "consensus-longranger/phased_variants.vcf.gz"

    if not os.path.isfile(longrangerBam) or not os.path.isfile(longrangerVCF):
        longrangerVCF, longrangerBam = tools.align_10x(consensusFa, queryReadsDir, outdir)
    else:
        print("Long Ranger files found, skipping step")

    #correct bug in longranger 2.2.2 where hets are sometimes given GT = 1|1 
        
    correctedVCF = re.sub(".gz", "", longrangerVCF)
    annotator.run_correction(longrangerVCF, correctedVCF)
    io.delete_file(longrangerVCF)
    
    longrangerVCF = tools.bgzip(correctedVCF)
    #longrangerVariants = [record for record in vcf.Reader(open(longrangerVCF, "rb"))]

    return longrangerBam

def reset_supplemental_alns(alignments, keepPrefix="hybrid"):
    
    # remove irrelevant SA 
    for alnmt in alignments:
        if alnmt.has_tag("SA"):
            sa = alnmt.get_tag("SA").split(";")
            filteredSA = []
            for a in sa:
                if not a.startswith(keepPrefix): continue
                filteredSA.append(a)
            
            if len(filteredSA) == 0:
                alnmt.set_tag("SA", None)
            else:
                alnmt.set_tag("SA", ";".join(filteredSA), "Z")
                
    alnDict = dict()
    for alnmt in alignments:
        if alnmt.qname not in alnDict:
            alnDict[alnmt.qname] = []
        alnDict[alnmt.qname].append(alnmt)

    for name in alnDict:
        hasMain = sum([(a.flag < 2048) for a in alnDict[name]]) > 0
        
        getqlen = lambda x: x.qlen
        if not hasMain:
            #print([a.qlen for a in alnDict[name]])
            longestAln = max(alnDict[name], key=getqlen)

            #set as non-supplemental
            if longestAln.flag >= 2048:
                longestAln.flag = longestAln.flag - 2048
                
    return alignments

def align_to_consensus_ref(consensusFa, refReads, rRegion, seqData, outdir, param):
    
    consensusAlignBam = outdir + "ref_aligned_consensus.pbmm2.bam"
    if os.path.isfile(consensusAlignBam):
        print("pbmm2 alignments found, skipping step")
        return consensusAlignBam

    #add flanking sequence to fasta file
    bufferLength = 10000
    
    sequenceDict = tools.fasta2dict(consensusFa, toUpper=True)
    seqName = list(sequenceDict.keys())[0]

    consensusRegion = SimpleRegion(seqName, 0, len(sequenceDict[seqName])-1)
    
    sequenceDict["l_flank"] = seqData[rRegion.chrom][max(0,rRegion.start-bufferLength) : rRegion.start]
    sequenceDict["r_flank"] = seqData[rRegion.chrom][rRegion.end+1 : rRegion.end + bufferLength]
    flankedFa = tools.dict2fasta(sequenceDict, outdir + "unpolished_plusflank")

    consensusFlankBam = tools.align_pacbio(flankedFa, refReads, outdir + "pacbio_aligned_consensus_plusflank")
    consensusAlignments = tools.samtools_fetch(consensusFlankBam, consensusRegion)
    consensusAlignments = reset_supplemental_alns(consensusAlignments, keepPrefix="hybrid")
    
    emptyBam = tools.samtools_write([], outdir + "temp", param.REF_ALIGNED_READS)
    bamHeader = tools.align_pacbio(consensusFa, emptyBam, outdir + "header")
    consensusAlignBam = tools.samtools_write(consensusAlignments, outdir + "ref_aligned_consensus.pbmm2", bamHeader)

    io.delete_file(consensusFlankBam)
    io.delete_file(emptyBam)
    io.delete_file(bamHeader)       

    return consensusAlignBam


    
