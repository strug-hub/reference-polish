import argparse

class Parameters:
    #default parameters established here
    def __init__(self):
        
        self.ID = "_id_missing_"
        self.OUTPUT_DIR   =   "./out"     #output directory
        self.REF_FA       =   None
        self.READS        =   None   
        self.REF_REGION   =   None   
        self.DOWNSAMPLE = 1

def get_parameters():
    """Parses command line arguments and sets default parameters.
    Returns a Parameters object."""

    p = Parameters()

    '''    
    p.ID = "CF002"
    p.REF_FA =  "/media/scott/Zapdos/reference/hg38.fa"
    p.READS = "/media/scott/Zapdos/slc9a3_polish/" + CFID + "_slc9a3.bam"
    p.OUTPUT_DIR = "/media/scott/Zapdos/slc9a3_polish/out/" + CFID

    #p.REF_REGION = "chr5:393462-696129"
    p.REF_REGION = "chr5:393462-677667"
    '''
    
    parser = argparse.ArgumentParser(description="Reference sequence region polish")

    # Positional
    parser.add_argument("refFa", metavar="reference", default=p.REF_FA, nargs="?",
                        help="Reference FASTA file")
    
    parser.add_argument("reads", metavar="reads", default=p.READS, nargs="?",
                    help="Pacbio reads")
    
    parser.add_argument("sampleId", metavar="id", default=p.ID, nargs="?",
             help="ID for sample")

    parser.add_argument("refRegion", metavar="region", default=p.REF_REGION, nargs="?",
                help="Region to polish, formatted as chrom:start-end")

    # Optional
    parser.add_argument("-o", "--outdir", type=str, default=p.OUTPUT_DIR,
                    help="Directory where output will be written." )
    
    parser.add_argument("-x", "--downsample", type=float, default=p.DOWNSAMPLE,
                    help="Randomly sample this propotion of the reads." )
    
    # Parse
    args = parser.parse_args()

    p.REF_FA = args.refFa
    p.REF_REGION = args.refRegion
    p.READS = args.reads
    p.ID = args.sampleId
    p.DOWNSAMPLE = args.downsample

    p.OUTPUT_DIR = args.outdir
    
    return p
