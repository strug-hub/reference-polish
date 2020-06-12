# reference-polish

Take pbmm2 aligned CLR reads and locally polish & phase the reference.

External Requirements
------

Tool | Location
--- | --- |
*samtools* | http://www.htslib.org/
*minimap2* | https://github.com/lh3/minimap2
*paftools.js* | https://github.com/lh3/minimap2/tree/master/misc
*pbmm2* | https://github.com/PacificBiosciences/pbmm2
*GenomicConsensus* | https://github.com/PacificBiosciences/gcpp
*longshot* | https://github.com/pjedge/longshot

Prior to running the script, please ensure all of these are installed and can be run. \
[external_tools.py](external_tools.py) assumes each tool can be run from the alias within this file. Ensure these alias exist on your system or replace the variables with the direct path to each tool.

For example: \
`MM2="/home/scott/bin/minimap2-2.17_x64-linux/minimap2"` \
`PAFTOOLS="k8 /home/scott/bin/minimap2-2.17_x64-linux/paftools.js"`

Running the script
------

`git clone https://github.com/ScottMastro/reference-polish.git` \
`cd reference-polish` \
`python main.py`

