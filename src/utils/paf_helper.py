import sys
sys.path.append("..")

import pandas as pd
from structures.region import SimpleRegion

def parse_paf(pafFile):
    """Loads file containing minimap2 alignments."""

    colnames = ["qid", "qlen", "qstart", "qend",
                "strand", 
                "rid", "rlen", "rstart", "rend", 
                "nmatch", "alnlen", "qual"]
    
    colorder = [c for c in colnames]
    
    typeMap = {"qlen": int, "qstart": int, "qend": int,
                    "rlen": int, "rstart": int, "rend": int,
                    "nmatch": int, "alnlen": int, "qual": int}
    
    lines = []
    with open(pafFile) as paf:
       line = paf.readline()
       while line:
           cols = line.split("\t")
           l = {c : x for c,x in zip(colnames, cols)}
           
           for i in range(len(colnames), len(cols)):
               col = cols[i]
               l[col[0:2]] = col[5:].rstrip("\n\r")
               if col[0:2] not in colorder: colorder.append(col[0:2])
               
           lines.append(l)
           line = paf.readline()
    
    
    df = pd.DataFrame(lines)
    df = df[colorder]
    
    if "dv" in df.columns:
        typeMap["dv"] = float    
    if "NM" in df.columns:
        typeMap["NM"] = int    


    df = df.astype(typeMap)

    return df

def paf_to_regions(pafFile, q=True):
    pafDf = parse_paf(pafFile)
    regions = []
    
    if q:
        for _,aln in pafDf.iterrows():
            regions.append(SimpleRegion(aln["qid"], aln["qstart"], aln["qend"]))
    else:
        for _,aln in pafDf.iterrows():
            regions.append(SimpleRegion(aln["rid"], aln["rstart"], aln["rend"]))

    return regions