This is a set of python scripts (and for some, the Conda environment and other tools requered to run them) to assemble, predict protein sequences in, and annotate the whole genome of a mollusc.  
The main pipeline is broken up into three parts:   
assembly.py which covers the long read de novo assembly from either Oxford Nanopore Technologies or Pac-Bio HiFi sequencing, an optional step of polishing the sequencing using short read genomic data, and clustering into chromosomes based on Hi-C data.  
After running assembly.py one must manually correct the Hi-C clustering until they recieve a satisfactory result, at which point the next script can be run.  
prediction.py covers the genome prediction using (!!!)  
annotation.py covers the functional genome annotation using InterProScan and BLAST+  
for more info, see the notebooks for the individual scripts.

A list of helper scripts is also included, these where used several times during the experimentation and might be useful for others looking to replicate or verify these results.  
These scripts, along with a brief description can be found in the helper/ directory
