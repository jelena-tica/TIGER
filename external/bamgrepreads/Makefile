SAMDIR=../samtools-0.1.17

BINDIR=.
INSTALLDIR=/usr/local/bin
CPP=g++
CC=gcc
CFLAGS=-O3 -Wall
COLOR_O=parsecolor.tab.o lex.parsecolor.o color.o
CHARSEQ=abstractcharsequence.o \
	genomicsequence.o \
	mutatedsequence.o \
	proteincharsequence.o \
	stringsequence.o
PARSESEG=segparser.tab.o lex.segparser.o segments.o


# edit: remove all targets except for bamgrepreads

TARGET_BINARIES=$(BINDIR)/bamgrepreads

all:$(TARGET_BINARIES)

$(BINDIR)/bamgrepreads: bamgrepreads.cpp
	$(CPP) -o $@ $^ $(CFLAGS) -I $(SAMDIR)  -L $(SAMDIR) -lbam   -lz




install:all
	cp $(TARGET_BINARIES) $(INSTALLDIR)

uninstall:
	rm -f $(subst $(BINDIR),$(INSTALLDIR),$(TARGET_BINARIES))

clean:
	rm -f $(TARGET_BINARIES) *.o a.out lex.segparser.c \
		segparser.tab.c segparser.tab.h segparser.output \
		lex.selectset.c parsecolor.tab.c lex.parsecolor.c  lex.yy.c  \
		selectset.tab.c parsecolor.tab.h selectset.tab.h \
		dbsnpxsl.cpp  genexsl.cpp  nuccorexsl.cpp  \
		pubmedxsl.cpp  taxonomyxsl.cpp selectset.output parsecolor.output

clean-all:clean uninstall
	

#test:$(CGI_BIN_DIR)/tabixserver.cgi
#test2: random.o tarball.o
#	g++ -Wall -g -o a.out vcfsim.cpp random.o tarball.o
test:$(BINDIR)/redondepth
	$<  -f /commun/data/pubdb/ucsc/hg19/chromosomes/hg19.fa \
		/commun/data/projects/NGS/201202macrogen/DataHD20120302/CD8262/CD8262_sorted.bam \
		/commun/data/projects/NGS/201202macrogen/DataHD20120302/CD5121/CD5121_sorted.bam \
		/commun/data/projects/NGS/201202macrogen/DataHD20120302/CD8260/CD8260_sorted.bam \
		/commun/data/projects/NGS/201202macrogen/DataHD20120302/CD8289/CD8289_sorted.bam \
		/commun/data/projects/NGS/201202macrogen/DataHD20120302/S0529/S0529_sorted.bam > /data/users/lindenb/jeter.redondepth.tsv

#test2:
#	 g++ -DTEST_THIS_CODE -Wall xstream.cpp `xml2-config --cflags --libs`
#	 ./a.out jeter.xml
