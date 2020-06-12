import os, re
import parameters
import utils.log as logger
import utils.fasta_handler as fasta

import polish.polish_region as polisher

def main():
    
    #--------------------------------------
    # READ IN DATA
    #--------------------------------------

    param = parameters.get_parameters_polish()
    if not os.path.exists(param.OUTPUT_DIR): os.mkdir(param.OUTPUT_DIR)
    
    logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)
    logger.Logger(level=param.VERBOSE, wait=param.WAIT)

    logger.Logger().out("Reading fasta...")
    
    if param.TARGET_CONTIG is None:
        seqData = fasta.read_fasta(param.FASTA)
    else:
        target = str(param.TARGET_CONTIG)
        seqData = {target : fasta.fasta_fetch(param.FASTA, target)}
        
    tigIds = list(seqData.keys())
    lengthData = {x : len(seqData[x]) for x in tigIds}
    tigIds.sort(key=lambda x: -lengthData[x])

    #logger.Logger().out("Reading alignment data...")
    #alignDict = io.parse_alignments(param.SUMMARY)

    #--------------------------------------
    # POLISH
    #--------------------------------------
    def clean_name(name): return re.sub(r'[\\/*?:"<>|]', "_", str(name))

    for tigId in tigIds:
        outdir = param.OUTPUT_DIR + "/" + clean_name(tigId) + "/"
        if not os.path.exists(outdir): os.mkdir(outdir)
        
        polisher.polish_contig(tigId, outdir, seqData, lengthData, param)
        
if __name__== "__main__":
  main()
  exit()
    
        
        
        
