import gzip

#the 16-base 10x barcode plus 7 additional bases
def trim_10x_barcode(fastq, trim=23):
    reader = gzip.open(fastq, "rt" if fastq.endswith("gz") else "r")
    trimmedFq = "/".join(fastq.split("/")[:-1]) +  "/trimmed" + str(trim) + "_" + fastq.split("/")[-1] 
    writer = gzip.open(trimmedFq, 'wt')

    line = reader.readline()

    while(line):        
        
        print(line)

        if not line.startswith("@") and not line.startswith("+"):
            line = line[trim:]
            print(line)
            
        writer.write(line)
        
        line = reader.readline()
        
    reader.close()
    writer.close()
    return trimmedFq
