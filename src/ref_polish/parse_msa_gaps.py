import argparse
import pandas as pd
import re
import os

def main():
    
    parser = argparse.ArgumentParser(description="Get repeat info from MSA")

    #Positional
    parser.add_argument("msa", metavar="msaFile", default="/media/scott/Rotom/hybrid2/slc9a3/out/mafft.align", nargs="?",
                        help="CLUSTAL format alignment file")
    parser.add_argument("repeats", metavar="repeatsFile", default="/media/scott/Rotom/hybrid2/slc9a3/out/repeats.txt", nargs="?", 
                    help="Tab-sep, 2 columns reference start/reference end positions; no header" )
    
    parser.add_argument("--out", type=str, default="/media/scott/Rotom/hybrid2/slc9a3/out/", nargs="?", 
                help="Output directory" )

    parser.add_argument("--refname", type=str, default="hg38", 
                    help="Reference name in MSA file")
    parser.add_argument("--start", type=int, default=0, 
                    help="Reference position for the first base in the MSA file")
    parser.add_argument("--buffer", type=int, default=10, 
                    help="# reference bases to include before and after each repeat")

    args = parser.parse_args()
    
    repeats = pd.read_csv(args.repeats, sep='\t', header=None)
    repeats = repeats.sort_values(by=0, ascending=True)
    
    sequence = dict()
    ref = args.refname
    
    with open(args.msa, "r") as reader:
        line = reader.readline()
        while line:
            
            line=re.sub("\s+", " ", line).strip()
            parts = line.split(" ")
            if len(parts) == 2 and parts[1].count("*") < 1:
                if parts[0] not in sequence:
                    sequence[parts[0]] = ""
                sequence[parts[0]] = sequence[parts[0]] + parts[1]
        
            line = reader.readline()
    reader.close()
    
    pos = 393462 #args.start
    posMap = dict()
    out = args.out + ("/" if args.out[-1] != "/" else "")
    outFasta = out + "fastas/"
                         
    if not os.path.exists(outFasta):
        os.mkdir(outFasta)

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