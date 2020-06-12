TOOL_DIR="/hpf/largeprojects/tcagstor/projects/cf_assembly_3g/workspace/mastros/hybrid/tools/"

CORES=16
MEM=121

#COMMON TOOLS
SAMTOOLS = "samtools"
BCFTOOLS = "bcftools"
BEDTOOLS = "bedtools"
BGZIP = "bgzip"
TABIX = "tabix"
PICARD="java -Xmx8G -jar " + TOOL_DIR + "picard.jar"

#https://anaconda.org/bioconda/nucdiff
#https://github.com/uio-cels/NucDiff
NUCDIFF="nucdiff"

#ILLUMINA TOOLS
BWA = "bwa"
PILON="java -Xmx16G -jar " + TOOL_DIR + "pilon-1.23.jar"

#PHASING TOOLS
WHATSHAP = TOOL_DIR + "whatshap"

HAP_PY = TOOL_DIR + "hap.py"
PRE_PY = TOOL_DIR + "pre.py"
RTG = TOOL_DIR + "rtg"

#GRAPH TOOLS
VG="/hpf/tools/centos7/vg/1.23.0/bin/vg"
NANOPLOT = TOOL_DIR + "NanoPlot"

#PACBIO TOOLS
MM2 = TOOL_DIR + "minimap2"
PBMM2 = TOOL_DIR + "pbmm2"
LONGSHOT = TOOL_DIR + "longshot"
PBSV = TOOL_DIR + "pbsv"
GENOMIC_CONSENSUS = TOOL_DIR + "gcpp"
CANU = "canu"

#10X TOOLS
LONGRANGER=TOOL_DIR + "longranger"
#required for longranger wgs
#version of GATK must be 3.3-3.8, or 4 except 3.6
GATK_10x_JAR=TOOL_DIR + "gatk-package-4.0.0.0-local.jar"

BAM2FQ_10X= TOOL_DIR + "bamtofastq-1.2.0"
