import subprocess
#import hpf_environment as env
import environment as env

def parse(command):
    return command.split()

def version(command, name, versionTag="--version"):
    
    if versionTag is not None and len(versionTag) > 0:
        command = command + " " + versionTag
    
    print("\nCHECKING: " + name + "\n" + command + "\n---------------")

    cmd = parse(command)
    
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    
    for x in stdout.decode('UTF-8').split("\n") + stderr.decode('UTF-8').split("\n"):
        if len(x) < 1: continue
        print(x)


if __name__== "__main__":
    print("cores = " + str(env.CORES))
    print("memory = " + str(env.MEM))

    print("ATTEMPTING TO IMPORT PYTHON PACKAGES...")
    import pysam
    import pyfaidx
    import pyfastx
    import subprocess
    import pybedtools
    import ssw
    import re
    import os
    from Bio import SeqIO
    import glob
    import shutil
    import gzip
    print("IMPORT SUCCESSFUL...")

    #COMMON TOOLS
    version(env.SAMTOOLS, "SAMTOOLS")
    version(env.BCFTOOLS, "BCFTOOLS")
    version(env.BEDTOOLS, "BEDTOOLS")
    version(env.BGZIP, "BGZIP")
    version(env.TABIX, "TABIX")
    #picard requires a tool name included to get version number
    version(env.PICARD, "PICARD", versionTag="ViewSam --version")

    #ILLUMINA TOOLS
    version(env.BWA, "BWA", versionTag=None)
    version(env.PILON, "PILON")
    
    
    #PHASING TOOLS
    version(env.WHATSHAP, "WHATSHAP")
    #version(env.HAP_PY, "HAP_PY")
    #version(env.PRE_PY, "PRE_PY")
    #print("PRE_PY = " + str(env.PRE_PY), "(exists = " + str(os.path.exists(env.PRE_PY)) + ")" )
    #version(env.RTG, "RTG")
    
    #GRAPH TOOLS
    version(env.VG, "VG", versionTag="version")
    version(env.NANOPLOT, "NANOPLOT")
    #version(env.GRAPHVIZ_DOT, "GRAPHVIZ_DOT", versionTag="-v")
    
    #PACBIO TOOLS
    version(env.MM2, "MM2")
    version(env.PBMM2, "PBMM2")
    version(env.LONGSHOT, "LONGSHOT")
    version(env.PBSV, "PBSV")
    version(env.GENOMIC_CONSENSUS, "GENOMIC_CONSENSUS")
    version(env.CANU, "CANU")
    
    #10X TOOLS
    version(env.LONGRANGER, "LONGRANGER", versionTag=None)
    print("Longranger properties:")
    print("GATK = " + str(env.GATK_10x_JAR), "(exists = " + str(os.path.exists(env.GATK_10x_JAR)) + ")" )
    version(env.BAM2FQ_10X, "BAM2FQ_10X", None)

    exit()
