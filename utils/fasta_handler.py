from Bio import SeqIO
import pyfaidx

import shutil
import re
import gzip
import csv
csv.field_size_limit(999999999)

rcDict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'C', 'c':'G', 'g':'C', 't':'A', 'N':'N', 'n':'N'} 
def reverse_complement(seq): return ''.join([rcDict[x] for x in seq[::-1]])

def index_fasta(faFile):
    """Index a FASTA file."""
    faidx = pyfaidx.Faidx(faFile)
    return faidx

def faidx_fetch(faFile, region):
    fa = pyfaidx.Faidx(faFile)
    return str(fa.fetch(region.chrom, region.start, region.end))

def fasta_fetch(faFile, fid):
    fa = pyfaidx.Fasta(faFile)
    return str(fa[fid])

def read_fasta(fasta, toUpper=False):
    """Loads fasta (or gzipped fasta) file into memory with SeqIO.
    Returns dictionary of fasta ID to sequence characters.
    """

    if(fasta[-2:] == "gz"):
        with gzip.open(fasta, "rt") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
    else:
        records = list(SeqIO.parse(fasta, "fasta"))

    def get_str(seq):
        return str(seq).upper() if toUpper else str(seq)
    
    faDict = dict(zip([str(r.id) for r in records], [get_str(r.seq) for r in records]))
    
    return faDict

#def write_fasta(fid, sequence, filePath):   
#    writer = open(filePath, "w+")
#    writer.write(">" + fid + "\n" + \
#             re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")

def write_fasta(faPath, fastaDict, toUpper=False, index=False):
    """Takes a directory path and dictionary of id->sequence.
    Writes a FASTA file from dictionary sequence.
    """

    if not faPath.endswith(".fa") and not faPath.endswith(".fasta"):
        faPath = faPath + ".fasta"
        
    writer = open(faPath, "w+")

    for fid in fastaDict:
        writer.write(">" + fid + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL) + "\n")

    writer.close()

    if index: index_fasta(faPath)
    return faPath

def get_fasta_len(faFile, fid=None):
    """Gets the sequence length of a FASTA sequence with id fid.
    Returns the first sequence length of fid is None.
    """
    faDict = read_fasta(faFile)
    
    if fid is None:
        fid = list(faDict.keys())[0]
    
    return len(faDict[fid])

def get_first_id(faFile):
    """Returns a string of the first id in the FASTA file.
    """

    faDict = read_fasta(faFile)
    fid = list(faDict.keys())[0]
    return fid

def get_fasta_seq(faFile, region=None, toUpper=False):
    """Returns a string of the sequence given by region.
    If region is None, returns the first FASTA sequence.
    """
    faDict = read_fasta(faFile, toUpper=toUpper)

    if region is None:
        fid = list(faDict.keys())[0]
        return faDict[fid]

    return faDict[region.chrom][region.start:region.end]


'''
def get_fasta_seq(faFile, region=None):   
    """Gets the sequence of a FASTA sequence with id fid.
    Returns the first sequence length of fid is None.
    """
    faDict = tools.fasta2dict(faFile)
    
    if region is None:
        fid = list(faDict.keys())[0]
        return faDict[fid]

    fid = region.chrom
    return faDict[fid][region.start:region.end]
'''

def rename_single_fasta(faFile, name, toUpper=False):
    seq = get_fasta_seq(faFile, toUpper=toUpper)
    faDict = {name : seq}
    temp = faFile + "__temporaryfa__.fasta"
    
    write_fasta(temp, faDict)
    shutil.move(temp, faFile)
    
    index_fasta(faFile)

    return faFile

def path_to_sequence(path, seqData):
    """
    Takes in a path, returns a list of sequences and a list of sources that represent
    the hybrid sequence when concatenated together.
    """

    sequence = []
    source = []
    
    def add_seq(startFork, endFork):
        
        if startFork.is_Nfork() or endFork.is_Nfork():
            sequence.append("N"*32)
            source.append("N"*32)
            return
        
        tigId = startFork.after_id()
        start = startFork.after_pos()
        end = endFork.before_pos()
        strand = startFork.after_strand()
        src = startFork.after_switch()        
        
        '''
        if(startFork.after_strand() != endFork.before_strand()):
            print(startFork)
            print(endFork)
            print("strand issue")
            input()
        '''
        
        seq = seqData[str(tigId)]
        
        if strand == -1:
            start = len(seq) - start
            end = len(seq) - end
            start, end = end, start
            
        segment = seq[start:end]
        if strand == -1:
            segment = reverse_complement(segment)

        sequence.append(segment)
        source.append(src*len(segment))
    
    startFork, endFork = None, path[0]
    
    for fork in path[1:]:
        startFork = endFork
        endFork = fork
        add_seq(startFork, endFork)
    
    return (sequence, source)

def write_hybrid(scaffolds, seqData, param):
    """
    Writes hybrid and source as FASTA.
    Also writes out source as a BED file.
    """
        
    hybridFasta = param.OUTPUT_DIR + "/hybrid_assembly.fasta"
    sourceFasta = param.OUTPUT_DIR + "/hybrid_source.fasta"
    sourceBed   = param.OUTPUT_DIR + "/hybrid_source.bed"
    
    sourceMap = {'r' : 'reference', 'q' : 'query', "N" : "NNN"}
    scoreMap = {'r' : '255', 'q' : '600', "N" : "999"}
    
    hf = open(hybridFasta, "w+")
    sf = open(sourceFasta, "w+")
    sb = open(sourceBed, "w+")

    for i, scaffold in enumerate(scaffolds):
        sequence, source = path_to_sequence(scaffold, seqData)
        tigId = scaffold.pid
        print("Writing " + tigId + "...")
        hf.write(">" + tigId + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL) + "\n")

        sf.write(">" + tigId + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(source), 0, re.DOTALL) + "\n")

        pos = 0
        for component in source:
            if len(component) > 0:
                x = component[0]
                sb.write("\t".join([tigId, str(pos), str(len(component)-1 + pos), \
                                    sourceMap[x], scoreMap[x], "."]) + "\n")
            pos = pos + len(component)
               
    hf.close()
    sf.close()
    sb.close()
    
def write_leftover(regions, seqData, param, minSize=10000, compressNs=10, spacerNs=500):
        
    fa = param.OUTPUT_DIR + "/hybrid_assembly_leftover.fasta"
    bed   = param.OUTPUT_DIR + "/hybrid_source_leftover.bed"
        
    f = open(fa, "w+")
    b = open(bed, "w+")

    fastaId="leftovers"
    f.write(">" + fastaId + "\n")

    pos=0
    sequences=[]
    for i, region in enumerate(regions):
        seq = seqData[region.chrom][region.start:region.end]
        split = seq.split('N')
        split = [x for x in split if len(x) > 0]
        if sum([len(x) for x in split]) < minSize:
            continue            

        seq = ("N"*compressNs).join(split)
        sequences.append(seq)
        
        b.write("\t".join([fastaId, str(pos), str(len(seq)-1 + pos), str(region)]) + "\n")
        pos += len(seq) + spacerNs
        
    sequence = ("N"*spacerNs).join(sequences)
    f.write(re.sub("(.{64})", "\\1\n", "".join(sequence), 0, re.DOTALL))
        
    f.close()
    b.close()
