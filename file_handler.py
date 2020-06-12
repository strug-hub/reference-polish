from pybedtools import BedTool
import shutil
import os
import re
import pickle as pkl
import pandas as pd
import pyfaidx

def make_dir(directory):
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

def file_exists(file):
    return os.path.isfile(file)

def delete_file(file, deleteDirTree=False):
    if os.path.isfile(file):
        os.remove(file)
    if os.path.isfile(file + ".bai"):
        os.remove(file + ".bai")
    if os.path.isfile(file + ".fai"):
        os.remove(file + ".fai")
    if os.path.isfile(file + ".tbi"):
        os.remove(file + ".tbi")

    if os.path.isdir(file):
        try:
            os.rmdir(file)
        except:
            if deleteDirTree: shutil.rmtree(file)

def rename_file(file, newFile):
    os.rename(file, newFile)
    return newFile


def faidx_fetch(faFile, region):
    fa = pyfaidx.Faidx(faFile)
    return str(fa.fetch(region.chrom, region.start, region.end))


def index_fasta(faFile):
    faidx = pyfaidx.Faidx(faFile)
    return faidx

def write_fasta(faPath, fastaDict, toUpper=False, index=False):

    if not faPath.endswith(".fa") and not faPath.endswith(".fasta"):
        faPath = faPath + ".fasta"
        
    writer = open(faPath, "w+")

    for fid in fastaDict:
        writer.write(">" + fid + "\n" + \
                 re.sub("(.{64})", "\\1\n", "".join(fastaDict[fid]), 0, re.DOTALL) + "\n")

    writer.close()

    if index: index_fasta(faPath)
    return faPath







import csv
csv.field_size_limit(999999999)

def copy_file(file, newFile):
    shutil.copyfile(file, newFile)
    return newFile

def move_file(file, newFile):
    shutil.move(file, newFile)
    return newFile






def parse_alignments(csv_file):
    """Loads file containing BLAST alignments."""

    colnames = ["qid", "qsize", "rid", "rlen",
                "chunks",  "rstart", "rend", 
                "qstart", "qend", "firstchunk", #todo: remove firstchunk col
                "nchunks", "alen", "pcid"]
    
    aligndf = pd.read_csv(csv_file, header=None, names=colnames, sep="\t", engine="python")
    alignDict = {str(tigId):rows for (tigId,rows) in tuple(aligndf.groupby('qid'))}

    return alignDict

def pickle(obj, path, overwrite=False):
    """Pickles an object.
    Optionally overwrite if file already exists at path."""

    exists = os.path.isfile(path)
    if overwrite or not exists:
         with open(path, 'wb') as handle:
            pkl.dump(obj, handle)
    
def unpickle(path):
    """Unpickles an object and return contents.
    Returns None if path does not exist"""

    exists = os.path.isfile(path)
    if exists:
        with open(path, 'rb') as handle:
            return pkl.load(handle)
    return None

def validate_ids(rids, qids):
    """Verifies that IDs between query and reference FASTA are unique.
    Returns True if IDs are ok, False if duplicate is found"""

    intersect = set(rids).intersection(qids) 
    if len(intersect) > 0:
        print("ERROR: IDs for input fasta files overlap.")
        print("\"" + str(intersect.pop()) + "\" is found in both fasta files.")
        print("Please ensure each ID is unique.")
        return False
    return True

def parse_bed(bedFile):
    """Loads BED file containing unitigs."""
    
    if bedFile is None:
        return None
    
    bed = BedTool(bedFile)        
    return bed
    
 
