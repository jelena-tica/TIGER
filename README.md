# TIGER
TIGER stands for Transduction Inference in GERmline genomes, 
and it performs mobile element mediated transduction detection.

The tool takes translocation prediction calls and mobile element 
insertion predictions and combines this information with the input 
from either a BWA or ELAND mapped BAM files to generate high-
confidence transduction callset.

Installing TIGER
----------------



Dependencies
------------

	* Blat >= v34
	* Bedtools (tested with version 2.17)
	* Perl >= 5.8.x

To install `blat`, visit the http://hgdownload.cse.ucsc.edu/admin/exe/, 
enter the folder correspoding to your machine type and download the 
blat pre-compiled binary (e.g. for Linux x86_64 machine, the following 
executable will be installed:
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat)

To install `bedtools` visit https://github.com/arq5x/bedtools2 and
download current the version (recommended). On a command line type:
  tar -zxvf BEDTools-<version>.tar.gz
  cd BEDTools-<version>
  make


