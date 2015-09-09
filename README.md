# TIGER
TIGER stands for **T**ransduction **I**nference in **GER**mline genomes, 
and it performs non-reference mobile element mediated transduction detection.

The tool takes mobile element insertion predictions and translocation 
prediction calls and  combines this information with the input 
from either a BWA or ELAND mapped BAMs files to generate a high-
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
ready to run. You best execute the script with its full pathname, e.g.

	/path/to/installation/TIGER/TIGER_v1.sh


Dependencies
------------

You need the following external tools installed and in your `PATH`:
* Blat >= v34
* Bedtools (tested with version 2.17)
* Perl >= 5.8.x

To install `blat`, visit http://hgdownload.cse.ucsc.edu/admin/exe/, 
enter the folder corresponding to your machine type and download the 
blat pre-compiled binary (e.g. for Linux x86_64 machine, the following 
executable will be installed:
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat)

To install `bedtools` visit https://github.com/arq5x/bedtools2 and
download the current version (recommended). On a command line type:
	
	tar -zxvf BEDTools-<version>.tar.gz
	cd BEDTools-<version>
	make

Run a test case
---------------
The folder `test` contains data from a small subset from an orangutan individual 
sequenced to 25x sequencing coverage in 2x100bp mode (Gokcumen et al., PNAS, 2013) that
you can use to test your installation. This data is aligned with Eland, the 
non-reference mobile element calls come from a modified version of TEA (Lee et al., Science, 2012)
and the translocation calls from Delly v0.0.11 - jumpy_v0.0.11 (Rausch et al., Bioinformatics, 2012).
You still need the reference genome FASTA file for orangutan, which can be downloaded from the [UCSC Genome Browser](http://hgdownload.soe.ucsc.edu/goldenPath/ponAbe2/bigZips/).

Then call TIGER like this (and replace `<ponAbe2>` with the Orangutan genome):

	./TIGER_v1.sh test/test.bam eland <ponAbe2.fa> test/test.l1 test/test.trans test/


Usage
-----
When you call TIGER without arguments a usage synopsis will be shown:

	TIGER_v1.sh bam aligner(bwa|eland) ref_genome L1_insertions translocations out_dir [min_match_len]

This input is required: 

* `bam`: Whole-Genome Sequencing BAM file
* `aligner`: State which aligner has been used to map the BAM file from above (either eland or bwa)
* `ref_genome`: FASTA file containing the reference genome (e.g. from UCSC)
* `L1_insertions`: L1 insertion calls (from your favorite MEI caller) in the 
  following TAB-separated format:

		chromosome_of_insertion	breakpoint1	breakpoint2	TSD

  where the `breakpoint1` and the `breakpoint2` are the start and end of a 
  Target-Site Duplication (TSD) on a `chromosome_of_insertion`.
* `translocations`: Translocation calls (from your favorite SV caller) in the 
  following TAB-separated format:

    	chromosome1	breakpoint1	chromosome2	breakpoint2

  where the translocation happened between two locations in the genome. 
  Both directions are taken into consideration (that the `chromosome1 breakpoint1` 
  translocates to `chromosome2 breakpoint2` and vice versa).
* `out_dir`: Folder to write results to
* `min_match_len` (optional): Specify the minimum length of BLAT alignments 
  to be considered. Defaults to 50. If you are using reads <=100bp, try lowering 
  this value.

The final transduction file will have the following format:
 
 	TargetChr       TargetStart     TargetEnd       SourceChr       SourceStart     SourceEnd       SumOfReads      AverageUniq     ReadsWith6A/Ts  TSD     SourceSize 

* `TargetChr, TargetStart, TargetEnd`: Target chromosome, target start, target end (locus of insertion)
* `SourceChr, SourceStart, SourceEnd`: Source chromosome, source start, source end (locus of origin)
* `SumOfReads`: Sum of reads indicating the translocation which confirms transduction
* `AverageUniq`: The average "uniqueness" of reads indicating transduction (mean of how many times the read is found in the genome)
* `ReadsWith6A/Ts`:  Number of reads containing at least six consecutive non-reference As/Ts
* `TSD`: Target Site Duplication (TSD) length
* `SourceSize`: Computationally assessed size of the predicted transduction



