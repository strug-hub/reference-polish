import argparse
import pandas as pd
import re
import os
import sys
import gzip

from Bio import SeqIO
import pyfaidx

def read_fasta(fasta, toUpper=False):
    """
    Loads fasta (or gzipped fasta) file into memory with SeqIO.
    Returns dictionary of id->sequence.
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


def main():
    
    parser = argparse.ArgumentParser(description="Get repeat info from MSA")

    #Positional
    parser.add_argument("msa", metavar="msaFile", nargs="?",
                        help="CLUSTAL format alignment file")
    parser.add_argument("repeats", metavar="repeatsFile", nargs="?", 
                    help="Tab-sep, 2 columns reference start/reference end positions; no header" )
    parser.add_argument("--out", type=str, default="./out/", nargs="?", 
                help="Output directory" )

    parser.add_argument("--refname", type=str, default="hg38", 
                    help="Reference name in MSA file")
    parser.add_argument("--start", type=int, default=0, 
                    help="Reference position for the first base in the MSA file (default=0)")
    parser.add_argument("--buffer", type=int, default=10, 
                    help="# reference bases to include before and after each repeat (default=10)")
    parser.add_argument("--samplefasta", type=str, default=None,
                         help="Extract and use sample names from this FASTA (because CLUSTAL format trims IDs)")
    
    args = parser.parse_args()
    
    if args.msa is None or args.repeats is None:
        parser.print_help()
        sys.exit()
    
    repeats = pd.read_csv(args.repeats, sep='\t', header=None)
    repeats = repeats.sort_values(by=0, ascending=True)
    
    sequence = dict()
    ref = args.refname
    
    with open(args.msa, "r") as reader:
        line = reader.readline()
        while line:
            if "*" in line:
                line = reader.readline()
                continue
                
            line=re.sub("\s+", " ", line).strip()
            parts = line.split(" ")
            if len(parts) == 2:
                if parts[0] not in sequence:
                    sequence[parts[0]] = ""
                sequence[parts[0]] = sequence[parts[0]] + parts[1]
        
            line = reader.readline()

    if args.samplefasta is not None:
        fasta = read_fasta(args.samplefasta)
        
        clustalKey=list(sequence.keys())
        newKey=list(fasta.keys())
        newSequence = {}
        for key in clustalKey:
            betterKey = [k for k in newKey if key.startswith(k)]
            if len(betterKey) != 1:
                print("ERROR: Could not match CLUSTAL id " + key + " uniquely to FASTA ID.")
                print("Possible options found: ", betterKey)
                exit()
            newSequence[betterKey[0]] = sequence[key]
        
        sequence=newSequence

    pos = args.start
    posMap = dict()
    out = args.out + ("/" if args.out[-1] != "/" else "")
    
    if not os.path.exists(out):
        os.mkdir(out)
    
    outFasta = out + "fastas/"
         
    if not os.path.exists(outFasta):
        os.mkdir(outFasta)

    if ref not in sequence:
        print("ERROR: Could not find reference ID " + ref + " in CLUSTAL file. Was it truncated?")
        print("IDs from CLUSTAL file: ", list(sequence.keys()))
        exit()
                
    for i,char in enumerate(sequence[ref]):
        if char != "-":
            posMap[pos] = i
            pos += 1
    
    samples = [sample for sample in sequence]
    samples.remove(ref) ; samples = sorted(samples)
    samples = [ref] + samples
    
    buffer = args.buffer
    
    countWriter= open(out + "counts.txt", "w")
    countWriter.write("start\tend\t" + "\t".join(samples) + "\n")
    
    for i,row in repeats.iterrows():
        s = row[0] - buffer
        e = row[1] + buffer
        name = str(s) + "_" + str(e)
        
        sIdx, eIdx = posMap[s], posMap[e]
        countData = [s, e]
        seqWriter = open(outFasta + name + ".fa", "w")
        seqGapWriter = open(outFasta + name + "_withgap.fa", "w")

        for sample in samples:
            seq = sequence[sample][sIdx:eIdx]
            seqNoGap = re.sub("-", "", seq)
            countData.append(len(seqNoGap))
            
            seqWriter.write(">" + sample + "_" + name + "\n" + \
                            re.sub("(.{64})", "\\1\n", "".join(seqNoGap), 0, re.DOTALL) + "\n")
            seqGapWriter.write(">" + sample + "_" + name + "\n" + \
                            re.sub("(.{64})", "\\1\n", "".join(seq), 0, re.DOTALL) + "\n")

        countWriter.write("\t".join([str(x) for x in countData]) + "\n")
        print(name, "MIN", min(countData[2:]), "MAX", max(countData[2:]))

        seqWriter.close()
        seqGapWriter.close()

    countWriter.close()


if __name__== "__main__":
  main() 
  print("done")
  #exit()
