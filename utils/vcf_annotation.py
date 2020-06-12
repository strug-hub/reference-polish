import argparse
import gzip
import pyfastx
import vcf

HEADERINFO="##vcf_correction=LongRanger 2.2.2 GT calls have been corrected\n" + \
            "##INFO=<ID=HetCorrected,Number=1,Type=Integer,Description=\"1 for variants that were corrected to het (LongRanger 2.2.2 bug)\">\n"

HEADERINFO_HP="##INFO=<ID=HOMP,Number=1,Type=Integer,Description=" + \
              "\"Homopolymer run in ref in bp (not including variant)\">\n" + \
              "##INFO=<ID=HOMP_DIFF,Number=1,Type=Integer,Description=" + \
              "\"Change in homopolymer count caused by variant\">\n"
              
CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE=0,1,2,3,4,5,6,7,8,9

'''
BARCODE_DICT = dict()

def get_barcodes(line):
    try:
        cols = line.rstrip().split("\t")
        sampDict = { f:s for f,s in zip(cols[FORMAT].split(":"), cols[SAMPLE].split(":")) }

        phaseSet = cols[CHROM] + "_" + sampDict["PS"]
        
        if phaseSet not in BARCODE_DICT:
            BARCODE_DICT[phaseSet] = {}
                
        if "|" not in sampDict["GT"]:
            return 
        
        A,B = [ int(x) for x in sampDict["GT"].split("|") ] 
        if A == B:
            return
        
        bxAllele = sampDict["BX"].split(",")
        
        for i,allele in enumerate(bxAllele):
            
            for barcode in allele.split(";"):
                bcInfo = barcode.split("-")
                bc, bcCount = bcInfo[0], bcInfo[1].count("_")

                hapCount = (bcCount if i == A else 0, bcCount if i == B else 0)
                if bc not in BARCODE_DICT[phaseSet]:
                    BARCODE_DICT[phaseSet][bc] = hapCount
                else:
                    prevCount= BARCODE_DICT[phaseSet][bc]
                    BARCODE_DICT[phaseSet][bc] = (hapCount[0] + prevCount[0], hapCount[1] + prevCount[1])       
    except:
        return 
    
    

def construct_barcode_dict(inputVCF):       
    if inputVCF.endswith(".gz"):
        reader = gzip.open(inputVCF, "rt")
    else:
        reader = open(inputVCF, "r")

    for line in reader:        
        if len(line) > 1 and line[0] != "#":
            get_barcodes(line)

    reader.close()
'''

def _annotate_hompolymers_line(line, pyfa, window=100):
    
    try:
        cols = line.split("\t")
        pos, ref, alt = int(cols[POS]), cols[REF], cols[ALT]

        #only annotate biallelic
        if "," in alt: return line

        #trim left
        while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
            pos += 1
            ref, alt = ref[1:], alt[1:]            
        #trim right
        while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]

        #only annotate indels
        if len(ref) == len(alt): return line
        
        #only annotate simple
        if len(ref) < 1:
            chars = set(alt)
        elif len(alt) < 1:
            chars = set(ref)
        else: return line
        
        if len(chars) != 1: return line
        nt = chars.pop()
        
        start = max(0, pos-window)
        region = pyfa.fetch(cols[CHROM], (start, pos+window))
        index = pos-start
        
        #print(ref,alt)
        #print(nt)
        #print(region[index-30:index], "|" , region[index:index+30])
        #print(region)

        count=0
        i = index-1
        while i >=0 and region[i]==nt:
            count += 1
            i -= 1
        
        i = index
        while i < len(region) and region[i]==nt:
            count += 1
            i += 1

        before, after = count, count + len(alt) - len(ref)

        #not a homopolymer
        if before == 0 or after == 0: return line
        #print(before,after)
        #print("HOMP="+str(before),"HOMPF="+str(after-before))
        
        return "\t".join(cols[:FORMAT]) + (";" if len(cols[INFO]) > 1 else "") + \
                "HOMP="+str(before) + ";HOMP_DIFF="+str(after-before) + "\t" + \
                 "\t".join(cols[FORMAT:])
    except:
        return line    
    
def annotate_hompolymers(inputVCF, outputVCF, refFa):

    headerWritten = False
    pyfa = pyfastx.Fasta(refFa)

    
    if inputVCF.endswith(".gz"):
        reader = gzip.open(inputVCF, "rt")
    else:
        reader = open(inputVCF, "r")

    if outputVCF.endswith(".gz"):
        writer = gzip.open(outputVCF, 'wt')
    else:
        writer = open(outputVCF, 'w')


    for line in reader:        

        if len(line) > 1:
            
            if not headerWritten:
                if line[0:2] == "##":
                    writer.write(line)
                elif line[0] == "#":
                    writer.write(HEADERINFO_HP)
                    writer.write(line)
                    headerWritten = True
                continue
                    
        annotatedLine = _annotate_hompolymers_line(line, pyfa)
        writer.write(annotatedLine)
        
    reader.close()
    writer.close()

    return outputVCF

def _fix_line(line):
    cols = line.split("\t")
    
    try:
        
        sampDict = { f:s for f,s in zip(cols[FORMAT].split(":"), cols[SAMPLE].split(":")) }
        #infos = cols[INFO].split(";")
        #infoDict = { v.split("=")[0]:v for v in infos }

        gtCount = int(sampDict["GT"][0]) + int(sampDict["GT"][-1])
        pl = [ int(x) for x in sampDict["PL"].split(",") ]
                
        if gtCount > 1 and pl[1] < pl[2]:
            print("Correcting " + cols[CHROM] + ":" + cols[POS], 
                  "\n\tGT=" + sampDict["GT"], "PL=" + str(pl))
            
            fixedLine = [ cols[CHROM], cols[POS], cols[ID], cols[REF], cols[ALT], cols[QUAL], cols[FILTER] ]
            fixedLine.append(cols[INFO] + ";" + "HetCorrected=1")
            sampDict["GT"] = "0/1"
            fixedLine.append(cols[FORMAT])
            samp = [ sampDict[f] for f in cols[FORMAT].split(":") ]
            fixedLine.append(":".join(samp))
        else:
            return line
        fixedLine = "\t".join(fixedLine)
        return fixedLine
    except:
        return line

def run_correction(inputVCF, outputVCF):
    
    #construct_barcode_dict(inputVCF)
    
    headerWritten = False
    
    if inputVCF.endswith(".gz"):
        reader = gzip.open(inputVCF, "rt")
    else:
        reader = open(inputVCF, "r")

    if outputVCF.endswith(".gz"):
        writer = gzip.open(outputVCF, 'wt')
    else:
        writer = open(outputVCF, 'w')


    for line in reader:        

        if len(line) > 1:
            
            if not headerWritten:
                if line[0:2] == "##":
                    writer.write(line)
                elif line[0] == "#":
                    writer.write(HEADERINFO)
                    writer.write(line)
                    headerWritten = True
                continue
                    
        correctedLine = _fix_line(line)
        writer.write(correctedLine)
        
    reader.close()
    writer.close()
    
    return outputVCF
    

def add_gt(vcfFile, newVCF):
    variants = [record for record in vcf.Reader(open(vcfFile, "r"))]
        
    reader = open(vcfFile, "r")
    writer = open(newVCF, "w+")

    while True:
        line = reader.readline()
        if line.startswith("##"):
            writer.write(line)
        else:
            break
    reader.close()
    
    writer.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")    
    header = ["#CHROM", "POS", "ID", "REF", "ALT",
               "QUAL",  "FILTER", "INFO", "FORMAT", "SAMPLE"]
    writer.write("\t".join(header) + "\n")
    
    for variant in sorted(variants, key=lambda x: x.POS):

        line = [variant.CHROM, variant.POS, ".", variant.REF, ",".join([str(x) for x in variant.ALT]), 
                variant.QUAL, "PASS", ",".join([str(key)+"="+str(variant.INFO[key]) for key in variant.INFO]),
                "GT", "0/1"]
        writer.write("\t".join([str(x) for x in line]) + "\n")

    writer.close()
    
    return newVCF

def combine_vcf_lines(vcf1, vcf2, newVCF):
    variants1 = [record for record in vcf.Reader(open(vcf1, "r"))]
    variants2 = [record for record in vcf.Reader(open(vcf2, "r"))]
    
    reader = open(vcf1, "r")
    writer = open(newVCF, "w+")

    while True:
        line = reader.readline()
        if line.startswith("##"):
            writer.write(line)
        else:
            break
    reader.close()
    
    #writer.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")    
    header = ["#CHROM", "POS", "ID", "REF", "ALT",
               "QUAL",  "FILTER", "INFO", "FORMAT", "SAMPLE"]
    writer.write("\t".join(header) + "\n")
    
    tups = [("0|1", v) for v in variants1]+[("1|0", v) for v in variants2]
    for gt,variant in sorted(tups, key=lambda x: x[1].POS):
                    
        line = [variant.CHROM, variant.POS, ".", variant.REF, ",".join([str(x) for x in variant.ALT]), 
                variant.QUAL, "PASS", ",".join([str(key)+"="+str(variant.INFO[key]) for key in variant.INFO]),
                "GT", gt]
        writer.write("\t".join([str(x) for x in line]) + "\n")

    writer.close()
    
    return newVCF

def filter_nonpass(vcfFile, newVCF):
    
    if vcfFile.endswith(".gz"):
        reader = gzip.open(vcfFile, "rt")
    else:
        reader = open(vcfFile, "r")

    variants = [record for record in vcf.Reader(reader)]
    filteredVariants = [r for r in variants if len(r.FILTER) < 1]
    reader.close()

    if vcfFile.endswith(".gz"):
        reader = gzip.open(vcfFile, "rt")
    else:
        reader = open(vcfFile, "r")
        
    writer = vcf.Writer(open(newVCF, "w+"), vcf.Reader(reader))
    
    for record in sorted(filteredVariants, key=lambda x: x.POS):
        writer.write_record(record)

    reader.close()
    writer.flush()
    writer.close()
    
    return newVCF


def main():
    
    parser = argparse.ArgumentParser(description="Longranger 2.2.2 VCF Correction Tool (single sample only)")
    parser.add_argument("vcf", metavar="vcf", default="phased_variants.vcf.gz", nargs="?", 
                    help="Longranger 2.2.2 \"phased_variants\" VCF" )
    args = parser.parse_args()


    inputVCF = args.vcf
    outputVCF = "/".join(inputVCF.split("/")[:-1]) + "/phased_variants.corrected.vcf"

    run_correction(inputVCF, outputVCF)


if __name__== "__main__":
  main()
  #exit()