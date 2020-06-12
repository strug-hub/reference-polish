import regions as rgn
import helper
import external_tools as tools
import vcf
import vcf_annotation as annotator
import gzip
import random
import math
import re
import os

Q,R,H,HG = "query", "ref", "hybrid", "hg38"

class AlignmentInfo:

    def _parse_cigar(self, aln):
        ops = "MIDNSHP=XB"
        cigars = aln.get_cigar_stats()
        return { op : cigars[0][i] for i,op in enumerate(ops) }       

    def _get_tag(self, tag, alns):
        try: 
            return [aln.get_tag(tag) for aln in alns]
        except:
            return None

    def __init__(self, alns):
        
        if not isinstance(alns, list):
            alns = [alns]
            
        alns = [aln for aln in alns if aln.reference_name is not None]

        self.alignments = len(alns)
        self.qlen = [aln.qlen for aln in alns]
        self.alen = [aln.alen for aln in alns]
        self.qsize = alns[0].infer_query_length() if len(alns) > 0 else 0
        self.strand = [-1 if aln.is_reverse else 1 for aln in alns]
        self.regions = [rgn.SimpleRegion(aln.reference_name, aln.reference_start, aln.reference_end) for aln in alns]                    

        self.cigars = [self._parse_cigar(aln) for aln in alns]
        self.compressedIndels = [aln.cigarstring.count("I") + aln.cigarstring.count("D") for aln in alns]
        self.concordance = self._get_tag("mc", alns)
        self.editDist = self._get_tag("NM", alns)

    def aligned_length(self):
        return sum(self.alen)
    def query_aligned_length(self):
        return sum(self.qlen)

    def indel_count(self):
        return sum([cig["I"] + cig["D"] for cig in self.cigars])
    def match_count(self):
        return sum([cig["M"] + cig["="] for cig in self.cigars])
    def mismatch_count(self):
        return sum([cig["X"] for cig in self.cigars])
    def clip_count(self):
        return sum([cig["S"] + cig["h"] for cig in self.cigars])

    def weighted_concordance(self):
        return sum([c*l for c,l in zip(self.concordance, self.qlen)])/self.query_aligned_length()

    def get_strand_count(self):
        plus, minus = 0, 0
        for s, alen in zip(self.strand, self.alen):
            plus += 0 if s < 0 else alen
            minus += 0 if s > 0 else alen
            
        return (plus, minus)
   
    def contiguous_region(self, misalignTolerance=0.95, sizeTolerance=1.25):        
        
        chromDict = {region.chrom: [] for region in self.regions}
        for region in self.regions: chromDict[region.chrom].append(region)
        totalSize = sum([len(region) for region in self.regions])
                
        #ignore small alignments to other chroms
        choice = None
        for chrom in chromDict: 
            percentage = sum([len(region) for region in chromDict[chrom]])/totalSize
            if  percentage > misalignTolerance:
                choice = chrom ; break
        
        if choice is None: return False
        chromRegion = chromDict[choice]
        chrom = choice

        start = min([region.start for region in chromRegion])
        end = max([region.end for region in chromRegion])
        rsize = end - start
        
        #ignore far-away small alignments
        while(True):

            if rsize > self.qsize*sizeTolerance:
                excludeStart = [region for region in chromRegion if region.start != start]
                excludeStartSize = max([region.end for region in excludeStart]) - min([region.start for region in excludeStart])
                excludeEnd = [region for region in chromRegion if region.end != end]
                excludeEndSize = max([region.end for region in excludeEnd]) - min([region.start for region in excludeEnd])
    
                #try removing start alignment
                if excludeStartSize < excludeEndSize:
                    s = sum([len(region) for region in excludeStart])
                    if s/totalSize > misalignTolerance:
                        chromRegion = excludeStart
                        start = min([region.start for region in chromRegion])
                        end = max([region.end for region in chromRegion])
                        rsize = end - start
                        if rsize >= self.qsize/sizeTolerance:
                            continue
                    
                #try removing end alignment
                s = sum([len(region) for region in excludeEnd])
                if s/totalSize > misalignTolerance:
                    chromRegion = excludeEnd
                    start = min([region.start for region in chromRegion])
                    end = max([region.end for region in chromRegion])
                    rsize = end - start
                    if rsize >= self.qsize/sizeTolerance:
                        continue

                return False
            
            else:
                break
            
        #print(rsize, self.qsize*tolerance, self.qsize/tolerance)
        if rsize <= self.qsize*sizeTolerance and rsize >= self.qsize/sizeTolerance:
            return rgn.SimpleRegion(chrom, start, end)
        return False
    
    #def __len__(self):
    #def __repr__(self):
    #def __str__(self):


def ref_genome_assessment(seqDict, outdir, param, requireContinuous, retryTolerance=None):
        
    fa = io.dict2fasta(outdir + "all_seq.fa", seqDict)
    bam = tools.align_bam(param.HG38_INDEX, fa,  outdir + "all_seq_to_hg38")
    alignments = tools.samtools_fetch(bam)
    
    summaries = dict()
    for x in seqDict:
        alns = [a for a in alignments if a.qname == x]
        summaries[x] = AlignmentInfo(alns)
    
    ok = True
    hgRegion = None
    minOverlap = 0.9
    for x in requireContinuous:
        region = summaries[x].contiguous_region()
        print("---------------")
        print(x, "alignments", summaries[x].regions)
        print(x, "continuous", region)
        if not region:
            if retryTolerance is not None:
                region = summaries[x].contiguous_region(sizeTolerance=retryTolerance)
            if not region:
                print(x + " did not contiguously align to reference genome.")
                print("alignments: ", summaries[x].regions)
    
                ok = False
                continue
            
        if hgRegion is None: hgRegion = region
        else:
            if region.chrom != hgRegion.chrom:
                print("Ambiguous chromosome mapping to reference genome.")
                print(hgRegion, "vs", region)
                ok = False
                break
            overlap = min(hgRegion.end, region.end) - max(hgRegion.start, region.start) 
            if overlap/len(region) < minOverlap or overlap/len(hgRegion) < minOverlap:
                print("Ambiguous position mapping to reference genome.")
                print(hgRegion, "vs", region)
                ok = False
                break
            hgRegion = rgn.SimpleRegion(
                    hgRegion.chrom,
                    min(hgRegion.start, region.start),
                    max(hgRegion.end, region.end))
    
    if not ok: hgRegion = None
    
    info = dict()
    
    info["hg_contiguous_alignment"] = hgRegion

    for x in seqDict:
        info["hg_align_length_" + x] = summaries[x].aligned_length()
        info["hg_align_count_" + x] = summaries[x].alignments
        info["hg_align_match_" + x] = summaries[x].match_count()
        info["hg_align_indel_" + x] = summaries[x].indel_count()
        info["hg_align_mismatch_" + x] = summaries[x].mismatch_count()
        
        region = summaries[x].contiguous_region()
        if not region and retryTolerance is not None:
            region = summaries[x].contiguous_region(sizeTolerance=retryTolerance)
        info["hg_contiguous_alignment_" + x] = region

    return (info, hgRegion)


def sequence_assessment(seqDict, param):
    
    def gap_count(seq):
        gapCompress = re.sub("[N]+", "N", seq)
        return gapCompress.count("N")
    def gc_count(seq):
        return (seq.count("G") + seq.count("C")) / (len(seq) - seq.count("N"))
    
    info = dict()
    for x in seqDict:
        
        info["length_" + x] = len(seqDict[x])
        info["ambiguous_bases_" + x] = seqDict[x].count("N")
        info["gaps_" + x] = gap_count(seqDict[x])
        info["gc_percent_" + x] = gc_count(seqDict[x])

    return info

def N50(alns):
    alns.sort(key=lambda x: x.alen, reverse=True)
    alenSum = sum([aln.alen for aln in alns])
    
    target = alenSum/2
    cumSum = 0
    index = None
    for i,aln in enumerate(alns):
        cumSum += aln.alen
        if cumSum > target:
            index = i
            break
        
    if index is None:
        return None
    return alns[index].alen
      
def make_read_dict(alns):
    readDict = {aln.qname : [] for aln in alns}
    for aln in alns:
        readDict[aln.qname].append(aln)
    return readDict
  
def pacbio_align_assessment(seqDict, rRegion, outdir, param):
    
    alignments = tools.samtools_fetch(param.REF_ALIGNED_READS, rRegion)
    
    def coverage(alignDict, length):
        return sum([sum([aln.alen for aln in alignDict[id]]) for id in alignDict]) / length
    
    readsBam = tools.samtools_write(alignments, outdir + "pacbio_alignments", param.REF_ALIGNED_READS)
    bams,alns,readDict = dict(), dict(), dict()

    for x in seqDict:
        xfa = io.dict2fasta(outdir + x + "_seq.fa", {x:seqDict[x]})
        bams[x] = tools.align_pacbio(xfa, readsBam, outdir + x + "_pacbio")
        alns[x] = tools.samtools_fetch(bams[x])
        readDict[x] = make_read_dict(alns[x])
        
    QUALITY_THRESHOLD = 80 #%
    SIZE_THRESHOLD = 30 #%
    EDGE_SIZE = 2000 #bp
    EDGE_THRESHOLD = 90 #%

    readQualityDict,sizePercentDict,isEdgeRead = dict(),dict(),dict()

    for x in readDict:
        for readId in readDict[x]:
            if readId not in readQualityDict:
                readQualityDict[readId] = []
                sizePercentDict[readId] = []
                isEdgeRead[readId] = []
                
            totalAlen = sum([aln.alen for aln in readDict[x][readId]])
            weightedQuality = sum([(aln.get_tag("mc")*aln.alen)/totalAlen for aln in readDict[x][readId]])
            percentAligned = sum([100*aln.qlen/aln.infer_read_length() for aln in readDict[x][readId]])

            start, end = 0, len(seqDict[x])
            startEdge, endEdge = EDGE_SIZE, max(0, end - EDGE_SIZE)


            startEdgeOverlap = sum([100*max(0, (min(startEdge, aln.reference_end) - max(start, aln.reference_start)))/EDGE_SIZE \
                               for aln in readDict[x][readId]])
            endEdgeOverlap = sum([100*max(0, (min(end, aln.reference_end) - max(endEdge, aln.reference_start)))/EDGE_SIZE \
                               for aln in readDict[x][readId]])
                        
            isEdgeOverlap = startEdgeOverlap > EDGE_THRESHOLD or endEdgeOverlap > EDGE_THRESHOLD

            readQualityDict[readId].append(weightedQuality)
            sizePercentDict[readId].append(percentAligned)
            isEdgeRead[readId].append(isEdgeOverlap)

    qualityFilterList = [readId for readId in readQualityDict if max(readQualityDict[readId]) < QUALITY_THRESHOLD ]
    sizeFilterList = [readId for readId in sizePercentDict if max(sizePercentDict[readId]) < SIZE_THRESHOLD ]
    edgeReadList = [readId for readId in isEdgeRead if sum(isEdgeRead[readId]) == len(isEdgeRead[readId]) ]

    filterSet = set(qualityFilterList).union(set(sizeFilterList) - set(edgeReadList))
    
    info = dict()
    for x in readDict:
        filteredReads = {readId : readDict[x][readId] for readId in readDict[x] if readId not in filterSet}
        filteredAlns = []
        for readId in filteredReads: filteredAlns.extend(filteredReads[readId])

        n50 = N50(filteredAlns)
        info["pacbio_align_n50_" + x] = 0 if n50 is None else n50
        info["pacbio_align_count_" + x] = len(filteredAlns)
        info["pacbio_read_count_" + x] = len(filteredReads)
        info["pacbio_read_count_" + x] = len(filteredReads)
        
        summaries = [AlignmentInfo(filteredReads[readId]) for readId in filteredReads]

        info["pacbio_align_bp_" + x] = sum([summary.query_aligned_length() for summary in summaries])
        info["pacbio_concordance_" + x] = 0 if info["pacbio_align_bp_" + x] == 0 else \
            sum([summary.weighted_concordance()*summary.query_aligned_length() for summary in summaries])/ \
                                            info["pacbio_align_bp_" + x]
            
        info["pacbio_match_" + x] = sum([summary.match_count() for summary in summaries])
        info["pacbio_indel_" + x] = sum([summary.indel_count() for summary in summaries])
        info["pacbio_mismatch_" + x] = sum([summary.mismatch_count() for summary in summaries])

    return info


def kmer_assessment(seqDict, outdir, param):
    
    k=21
    
    def assess(compress):
        kmers = dict()
        for x in seqDict:
            kmers[x] = kt.KmerSet(seqDict[x], k, compressed=compress)
        
        seqList = [x for x in seqDict]
        for i,x in enumerate(seqList):
            for y in seqList[i+1:]:
                pre = "compressed_" if compress else ""
                info[pre + "kmer_similarity_" + x + "_" + y] = kmers[x].jaccard_similarity(kmers[y])
    
        if len(seqList) == 3:
            kmerBreakdown = kt.draw_pie3(seqList[0], kmers[seqList[0]],
                                         seqList[1], kmers[seqList[1]],
                                         seqList[2], kmers[seqList[2]],
                                         outdir + pre + "kmer_breakdown")
        elif len(seqList) == 2:
            kmerBreakdown = kt.draw_pie2(seqList[0], kmers[seqList[0]],
                             seqList[1], kmers[seqList[1]],
                             outdir + pre + "kmer_breakdown")
            
        return kmerBreakdown
    
    
    info = dict()
    info.update(assess(False))
    info.update(assess(True))

    return info

def write_info(info, outdir, outputFileName, transpose=True):
    
    order = [key for key in info] ; order.sort()
    values = [str(info[x]) for x in order]
    
    out = open(outdir + outputFileName, 'w')
    out.write("\t".join(order) + "\n")
    out.write("\t".join(values))
    out.close()
    
    if transpose:
        out = open(outdir + "transpose_" + outputFileName, 'w')    
        out.write("\n".join([o + ":\t" + v for o,v in zip(order,values)]))    
        out.close()
        
def get_outdir(prefix, param, fork1, fork2):
    
    outdir = param.OUTPUT_DIR + "/" + prefix + "_" + fork1.qid + "_" + \
                              str(fork1.qpos) + "_" + \
                              str(fork2.qpos) + "/"
                              
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    return outdir
    
def get_regions(fork1, fork2, lengthData):
    
    qid, rid = fork1.qid, fork1.rid
    if qid != fork2.qid: 
        qRegion = None
    else:
        qp = [fork1.get_pos_norm(lengthData, q=True), 
              fork2.get_pos_norm(lengthData, q=True)]
        qRegion = rgn.SimpleRegion(qid, min(qp[0], qp[1]), max(qp[0], qp[1]))
        
    if rid != fork2.rid:
        rRegion = None
    else:
        rp = [fork1.get_pos_norm(lengthData, q=False), 
              fork2.get_pos_norm(lengthData, q=False)]
        rRegion = rgn.SimpleRegion(rid, min(rp[0], rp[1]), max(rp[0], rp[1]))

    return (qRegion, rRegion)

def analyze_sequence(blockPath, seqData, lengthData, param, skipIfDone=True):

    seqFa = ["consensus.fa"]
    
    

def analyze_block(blockPath, seqData, lengthData, param, skipIfDone=True):
    
    
    
    
    outputFileName = "info.txt"
    info = dict()

    outdir = get_outdir("block", param, blockPath[0], blockPath[-1])
    
    print("Dumping info to " + outdir)
    if skipIfDone and os.path.isfile(outdir + outputFileName):
        print("Analysis has already been done, skipping.")
        return
    
    qRegion, rRegion = get_regions(blockPath[0], blockPath[-1], lengthData) 
    qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
    rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]
    hSeq, source = io.path_to_sequence(blockPath, seqData)
    hSeq, source = "".join(hSeq), "".join(source)
    info["hybrid_r_composition"] = source.lower().count("r")/len(source)

    print("LENGTHS:", len(qSeq), len(rSeq), len(hSeq))
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    seqSummary = sequence_assessment(seqDict, param)
    info.update(seqSummary)
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    refGenomeSummary, hgRegion = ref_genome_assessment(seqDict, outdir, param, requireContinuous=[Q,R])

    info.update(refGenomeSummary)
    
    if hgRegion is None:
        return
        print(" WAIT! ")
        input()   
    else: print(hgRegion)

    hgSeq = tools.faidx_fetch(param.HG38, hgRegion) if hgRegion is not None else ""
    print("hybrid vs reference length", len(hSeq), len(hgSeq))

    seqDict = { HG: hgSeq }
    hgSeqSummary = sequence_assessment(seqDict, param)
    info.update(hgSeqSummary)

    seqDict = { Q: qSeq, R: rSeq, HG: hgSeq }
    kmerSummary = kmer_assessment(seqDict, outdir, param)
    info.update(kmerSummary)

    seqDict = { Q: qSeq, R: rSeq, H: hSeq, HG: hgSeq }
    pbSummary = pacbio_align_assessment(seqDict, rRegion, outdir, param)
    info.update(pbSummary)
    
    write_info(info, outdir, outputFileName)

    
'''
def analyze_block(blockPath, seqData, lengthData, param, skipIfDone=True):
    
    outputFileName = "info.txt"
    info = dict()

    outdir = get_outdir("block", param, blockPath[0], blockPath[-1])
    
    print("Dumping info to " + outdir)
    if skipIfDone and os.path.isfile(outdir + outputFileName):
        print("Analysis has already been done, skipping.")
        return
    
    qRegion, rRegion = get_regions(blockPath[0], blockPath[-1], lengthData) 
    qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
    rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]
    hSeq, source = io.path_to_sequence(blockPath, seqData)
    hSeq, source = "".join(hSeq), "".join(source)
    info["hybrid_r_composition"] = source.lower().count("r")/len(source)

    print("LENGTHS:", len(qSeq), len(rSeq), len(hSeq))
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    seqSummary = sequence_assessment(seqDict, param)
    info.update(seqSummary)
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    refGenomeSummary, hgRegion = ref_genome_assessment(seqDict, outdir, param, requireContinuous=[Q,R])

    info.update(refGenomeSummary)
    
    if hgRegion is None:
        return
        print(" WAIT! ")
        input()   
    else: print(hgRegion)

    hgSeq = tools.faidx_fetch(param.HG38, hgRegion) if hgRegion is not None else ""
    print("hybrid vs reference length", len(hSeq), len(hgSeq))

    seqDict = { HG: hgSeq }
    hgSeqSummary = sequence_assessment(seqDict, param)
    info.update(hgSeqSummary)

    seqDict = { Q: qSeq, R: rSeq, HG: hgSeq }
    kmerSummary = kmer_assessment(seqDict, outdir, param)
    info.update(kmerSummary)

    seqDict = { Q: qSeq, R: rSeq, H: hSeq, HG: hgSeq }
    pbSummary = pacbio_align_assessment(seqDict, rRegion, outdir, param)
    info.update(pbSummary)
    
    write_info(info, outdir, outputFileName)

    
def analyze_inter_block(prevBlockPath, blockPath, seqData, lengthData, param, skipIfDone=True):

    outputFileName = "info.txt"
    info = dict()

    outdir = get_outdir("inter", param, prevBlockPath[-1], blockPath[0])

    print("Dumping info to " + outdir)
    if skipIfDone and os.path.isfile(outdir + outputFileName):
        print("Analysis has already been done, skipping.")
        return

    qid, rid = prevBlockPath[-1].qid, prevBlockPath[-1].rid
    if qid != blockPath[0].qid or rid != blockPath[0].rid:
        print("Blockpath IDs don't match...")
    
    qRegion, rRegion = get_regions(prevBlockPath[-1], blockPath[0], lengthData) 
    qSeq = seqData[qRegion.chrom][qRegion.start:qRegion.end]
    rSeq = seqData[rRegion.chrom][rRegion.start:rRegion.end]
    hSeq, source = rSeq, "r"*len(rSeq)
    info["hybrid_r_composition"] = source.lower().count("r")/len(source)

    print("LENGTHS:", len(qSeq), len(rSeq), len(hSeq))
    
    seqDict = { Q: qSeq, R: rSeq, H: hSeq }
    seqSummary = sequence_assessment(seqDict, param)
    info.update(seqSummary)
    
    seqDict = { Q: qSeq, R: rSeq }
    refGenomeSummary, hgRegion = ref_genome_assessment(seqDict, outdir, param, requireContinuous=[R], retryTolerance=5)
    info.update(refGenomeSummary)
    
    if hgRegion is None:
        return
        print(" WAIT! ")
        input()   
    else: print(hgRegion)

    hgSeq = tools.faidx_fetch(param.HG38, hgRegion) if hgRegion is not None else ""
    print("hybrid vs reference length", len(hSeq), len(hgSeq))
    
    seqDict = { HG: hgSeq }
    hgSeqSummary = sequence_assessment(seqDict, param)
    info.update(hgSeqSummary)

    seqDict = { R: rSeq, HG: hgSeq }
    kmerSummary = kmer_assessment(seqDict, outdir, param)
    info.update(kmerSummary)

    seqDict = { Q: qSeq, R: rSeq, HG: hgSeq }
    pbSummary = pacbio_align_assessment(seqDict, rRegion, outdir, param)
    info.update(pbSummary)

    write_info(info, outdir, outputFileName)
'''