SAMDIR=samtools-0.1.17
BAMGREPDIR=bamgrepreads

.PHONY: all clean

all:
	(cd ${SAMDIR}; make)
	(cd ${BAMGREPDIR};make)

clean:
	(cd ${SAMDIR}; make clean)
	(cd ${BAMGREPDIR}; make clean)

