#COMMON TOOLS
SAMTOOLS = "samtools"
BCFTOOLS = "bcftools"
BEDTOOLS = "bedtools"
BGZIP = "bgzip"
TABIX = "tabix"
PICARD="java -Xmx8G -jar /home/scott/bin/picard.jar"

CORES=8
MEM=16

#https://anaconda.org/bioconda/nucdiff
#https://github.com/uio-cels/NucDiff
NUCDIFF="nucdiff"

#ILLUMINA TOOLS
BWA = "bwa"
PILON="java -Xmx16G -jar /home/scott/bin/pilon-1.23.jar"

#PHASING TOOLS
WHATSHAP = "whatshap"
HAP_PY="/home/scott/bin/hap.py-install/bin/hap.py"
PRE_PY="/home/scott/bin/hap.py-install/bin/pre.py"
RTG = "/home/scott/bin/rtg-tools-3.10.1/rtg"

#GRAPH TOOLS
VG="/home/scott/bin/vg"
NANOPLOT = "NanoPlot"
GRAPHVIZ_DOT="dot"

#PACBIO TOOLS
MM2 = "/home/scott/bin/minimap2-2.17_x64-linux/minimap2"
PAFTOOLS = "k8 /home/scott/bin/minimap2-2.17_x64-linux/paftools.js"

PBMM2 = "pbmm2"
LONGSHOT = "longshot"
PBSV = "pbsv"
GENOMIC_CONSENSUS = "gcpp"
CANU = "/home/scott/bin/canu-1.9/Linux-amd64/bin/canu"

#10X TOOLS
LONGRANGER="/home/scott/bin/longranger-2.2.2/longranger"
#required for longranger wgs
#version of GATK must be 3.3-3.8, or 4 except 3.6
GATK_10x_JAR="/home/scott/bin/gatk-4.0.0.0/gatk-package-4.0.0.0-local.jar"

BAM2FQ_10X="/home/scott/bin/bamtofastq-1.2.0"
