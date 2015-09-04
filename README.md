# TIGER
TIGER stands for **T**ransduction **I**nference in **GER**mline genomes, 
and it performs mobile element mediated transduction detection.

The tool takes translocation prediction calls and mobile element 
insertion predictions and combines this information with the input 
from either a BWA or ELAND mapped BAM files to generate high-
confidence transduction callset.

Installing TIGER
----------------
TIGER consists of a single bash script that mostly uses standard Unix 
commands. However, it relies on certain external tools that have to be 
pre-installed (see Dependencies) or are packaged into this repository
for your convenience and to avoid problems with conflicting versions.
Included are `bamgrepreads`, a tool from [Pierre Lindenbaum](http://code.google.com/p/variationtoolkit),
and a version of [samtools](http://sourceforge.net/projects/samtools/files/samtools/0.1.17/). 

Obtain this repository to your computer by either cloning it

	git clone https://github.com/jelena-tica/TIGER.git

or downloading it as a zip file. Uncompress the zip in the latter case.

Enter the `external` directory

	cd TIGER/external

and simply type

	make

This will build samtools and bamgrepreads. Afterwards TIGER should be 
ready to go. You best the script with its full pathname, e.g.

	/path/to/installation/TIGER/TIGER_v1.sh


Dependencies
------------

You need the following external tools installed and in your `PATH`:
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




