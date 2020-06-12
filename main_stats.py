import utils.log as logger
import utils.parameters as parameters
import utils.file_handler as io
import utils.fasta_handler as fasta

import stitch.stitcher as stitcher
import weld.welder as welder
import scaffold.scaffolder as scaffolder

def main_stats():

    #--------------------------------------
    # READ IN DATA
    #--------------------------------------

    param = parameters.get_parameters()

    logger.FileLogger(clean=True, outdir=param.OUTPUT_DIR)
    logger.Logger(level=param.VERBOSE, wait=param.WAIT)

    logger.Logger().out("Reading Canu fasta...")
    refData = fasta.read_fasta(param.REF_FA)
    logger.Logger().out("Reading Supernova fasta...")
    queryData = fasta.read_fasta(param.QUERY_FA)
    
    seqData = dict()
    seqData.update(refData) ; seqData.update(queryData)
    lengthData = {x : len(seqData[x]) for x in seqData.keys()}
    
    logger.Logger().out("Reading alignment data...")
    alignDict = io.parse_alignments(param.SUMMARY)

    rids, qids = list(refData.keys()), list(queryData.keys())
    #if not io.validate_ids(rids, qids): return
    rids.sort(key=lambda x: -lengthData[x])
    qids.sort(key=lambda x: -lengthData[x])

    confBed = io.parse_confident_regions(param.REF_BED, rids, param)   
    
    #--------------------------------------
    # STITCH BLOCKS
    #--------------------------------------

    contigs = io.unpickle(param.OUTPUT_DIR + "/contigs.pickle")

    stitcher.stitch_stats()        
    contigs = []
    for tigId in qids:
        logger.Logger().out("Stitching contig: " + tigId, 1)
        contig = stitcher.stitch(tigId, alignDict, lengthData, param)
        if contig is not None: contigs.append(contig)

    logger.FileLogger().flush_all()
    io.pickle(contigs, param.OUTPUT_DIR + "/contigs.pickle")


    #todo: output information about BLOCK ORDER! (missassembly detection)


    #--------------------------------------
    # WELD BLOCKS
    #--------------------------------------
    
    paths = io.unpickle(param.OUTPUT_DIR + "/paths.pickle")
    if paths is not None:    
        logger.Logger().out("Loading welded data...")
    else:
        logger.Logger().out("Welding contigs...")

        paths = []
        for contig in contigs:
            path = welder.weld(contig, seqData, lengthData, param)
            if path is not None: paths.append(path)

        logger.FileLogger().flush_all()
        io.pickle(paths, param.OUTPUT_DIR + "/paths.pickle")        
        
    #--------------------------------------
    # SCAFFOLD BLOCKS
    #--------------------------------------

    rawScaffolds = io.unpickle(param.OUTPUT_DIR + "/scaffolds_raw.pickle")
    if rawScaffolds is not None:    
        logger.Logger().out("Loading scaffolded data...")
    else:
        logger.Logger().out("Scaffolding contigs...")

        ridSet = set()
        for path in paths:
            for fork in path: ridSet.add(fork.rid)
        
        for tigId in rids:
            if tigId not in ridSet: continue
            logger.Logger().out("Scaffolding with contig: " + tigId, 1)
            paths, exclude = scaffolder.scaffold(paths, tigId, confBed, lengthData, param)

        rawScaffolds = paths
        logger.FileLogger().flush_all()
        io.pickle(rawScaffolds, param.OUTPUT_DIR + "/scaffolds_raw.pickle")        

    #--------------------------------------
    # SALVAGE LEFTOVERS
    #--------------------------------------

    scaffolds = io.unpickle(param.OUTPUT_DIR + "/scaffolds.pickle")
    if scaffolds is not None:    
        logger.Logger().out("Loading salvaged data...")
    else:
        logger.Logger().out("Salvaging contigs...")

        scaffolds, unusedDict = scaffolder.salvage(rawScaffolds, qids, rids, lengthData, minSize=10000)   
        io.pickle(scaffolds, param.OUTPUT_DIR + "/scaffolds.pickle")

    #--------------------------------------
    # OUTPUT
    #--------------------------------------

    logger.Logger().out("Writing sequence...")

    fasta.write_hybrid(scaffolds, seqData, param)
    
    leftoverRegions = []
    for key in unusedDict:
        leftoverRegions.extend(unusedDict[key])
    fasta.write_leftover(leftoverRegions, seqData, param,
                      minSize=10000, compressNs=10, spacerNs=500)
    
    logger.FileLogger().flush_all()


if __name__== "__main__":
  main_stats()
  exit()