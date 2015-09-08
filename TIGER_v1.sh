#!/bin/bash

if [ $# -lt 6 ]
then
	echo -e "Usage:\n`basename $0` bam aligner(bwa|eland) ref_genome L1_insertions translocations out_dir [min_match_len]\n"
	exit 1
fi

MINMATCHLEN=50
if [ $# -gt 6 ] && [[ $7 =~ ^-?[0-9]+$ ]]
then
	MINMATCHLEN=$7
fi	


#BINDIR=$(readlink -f $0)
#BINDIR=${BINDIR%/*}
BINDIR=$(dirname $0)
echo "The script is in $BINDIR"

# bundle dependencies
SAMTOOLS="${BINDIR}/external/samtools-0.1.17/samtools"
REFORMATBLAT="${BINDIR}/external/reformatingBlat.pl"
BAMGREPREADS="${BINDIR}/external/bamgrepreads/bamgrepreads"
SP2TAB="${BINDIR}/external/sp2tab.py"
TAB2SP="${BINDIR}/external/tab2sp.py"

# tool dependencies
BEDTOOLS="bedtools"
BLAT="blat"

# I/O
BAM="$1"
ALIGNER="$2"
REFGENOME="$3"
REFNAME=$(basename $REFGENOME)
L1INS=$4
L1INSNAME=$(basename $L1INS)
TRANS=$5
TRANSNAME=$(basename $TRANS)
OUTDIR="${6%/}"
OUTPUT="$OUTDIR/transductions.bed"
SAMPLE=$(basename $BAM .bam)
TMPDIR="${6%/}/tmp"

 
echo "--> Configuration:"
echo "    This script is hosted in: $BINDIR"
if [ ! -f "$BAM" ]; then
	echo "    Error: Input BAM file not found: $BAM"
	exit 1
else
	echo "    Using BAM file:           $BAM"
fi
if [ "$ALIGNER" == "bwa" ] || [ "$ALIGNER" == "eland" ]; then
	echo "    Specified aligner:        $ALIGNER"
else
	echo "    Error: Aligner must be either 'eland' or 'bwa'"
	exit 1
fi
if [ ! -f "$REFGENOME" ]; then
	echo "    Error: Input ref genome file not found: $REFGENOME"
	exit 1
else
	echo "    Using reference genome:   $REFGENOME"
fi
if [ ! -f "$L1INS" ]; then
	echo "    Error: Input L1 insertion file not found: $L1INS"
	exit 1
else
	echo "    Using L1 insertion file:  $L1INS"
fi
if [ ! -f "$TRANS" ]; then
	echo "    Error: Input translocation file not found: $TRANS"
	exit 1
else
	echo "    Using translocation file: $TRANS"
fi
if [ ! -d "$OUTDIR" ]; then
	echo "    Error: The output directory does not exist: $OUTDIR"
	exit 1
else
	echo "    The output dir is:        $6"
fi
echo "    The output file is:       $OUTPUT"
mkdir $TMPDIR 2> /dev/null
if [ ! -d "$TMPDIR" ]; then
	echo "    Error: Could not create TMP dir: $TMPDIR"
	exit 1
else
	echo "    I created a tmp dir at:   $TMPDIR"
fi



echo "--> Preparing files for TIGER ..."
#1. take all L1 calls and reformat start-end (+/-500bp)

cat ${L1INS} | awk '{print $1"\t"$2-500"\t"$3+500}' > ${TMPDIR}/${L1INSNAME}.bed 

#2. take translocation calls and reformat start-end (+1bp)

cat ${TRANS} | awk '{print $1"\t"$2"\t"$2+1}' > ${TMPDIR}/${TRANSNAME}.first.bed
cat ${TRANS} | awk '{print $3"\t"$4"\t"$4+1}' > ${TMPDIR}/${TRANSNAME}.second.bed

echo "--> Initiating overlap of L1 calls and translocation calls ..."
#3. overlap L1 calls / translocation calls - to filter out L1 calls that are not overlapping translocations

#cd ${TMPDIR}
for i in first second
do
	${BEDTOOLS} intersect -a ${TMPDIR}/${L1INSNAME}.bed -b ${TMPDIR}/${TRANSNAME}.${i}.bed -wa -wb \
	> ${TMPDIR}/${SAMPLE}.${i}Coord 
	cut -f 1,2,3 ${TMPDIR}/${SAMPLE}.${i}Coord  \
	| sort \
	| uniq > ${TMPDIR}/${SAMPLE}.${i}Coord.regions 

done


echo "--> Extracting all translocation and single-anchored reads and mates from a BAM file ..."
#4. pulling out reads in each region 

for i in first second
do
	${BEDTOOLS} intersect -abam ${BAM} -b  ${TMPDIR}/${SAMPLE}.${i}Coord.regions >  ${TMPDIR}/${SAMPLE}.${i}ExtractedRegions.bam 
	${SAMTOOLS} view  ${TMPDIR}/${SAMPLE}.${i}ExtractedRegions.bam | cut -f 1 >  ${TMPDIR}/${SAMPLE}.${i}ExtractedRegions.bam.regionID 
done


#5. getting the mates
#http://code.google.com/p/variationtoolkit/

for i in first second
do
	${BAMGREPREADS} -R ${TMPDIR}/${SAMPLE}.${i}ExtractedRegions.bam.regionID  ${BAM} > ${TMPDIR}/${SAMPLE}.${i}bothMates 


	if [ "$ALIGNER" = "bwa" ]
	then
		cat ${TMPDIR}/${SAMPLE}.${i}bothMates  \
		| awk '{if ($7!="=") {print $0}}' \
		| sort -k7,8 \
		> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLmate 

		${SAMTOOLS} view -f 0x0080 ${TMPDIR}/${SAMPLE}.${i}ExtractedRegions.bam \
		| awk '{if ($4==$8) {print $0}}' \
		| sort -k7,8 \
		> ${TMPDIR}/${SAMPLE}.${i}bothMates.SAmate 

		cat ${TMPDIR}/${SAMPLE}.${i}bothMates.TLmate ${TMPDIR}/${SAMPLE}.${i}bothMates.SAmate > ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate 
	elif [ "$ALIGNER" = "eland" ]
        then
        cat ${TMPDIR}/${SAMPLE}.${i}bothMates  \
		| awk '{if (($7!="=")&&($7!="*")) print $0}' \
		| sort -k7,8 \
		> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate 
	fi
	awk '{print $1 "_" $2 "\t" $0}' ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate  | sort -k 1,1 > ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted
done

echo "--> Preparing for realignment step with BLAT ..."
for i in first second
do
	${SAMTOOLS} view  ${TMPDIR}/${SAMPLE}.${i}ExtractedRegions.bam \
	| awk '{print $1 "_" $2 "\t" $0}' \
	| sort -k1,1 \
	| join -v2 -1 1 -2 1 -  ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted \
	| python ${SP2TAB} \
	| cut -f2- \
	>  ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat
done


for i in first second
do
	cat ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat | cut -f 1,2,3,4,7,8,10 | awk '{print ">"$1"_"$2"_"$3"_"$4"_"$5"_"$6"\n"$7}' > ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa 
	cat ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat | cut -f 1,2,3,4,7,8,10 | awk '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"\t"$1"\t"$2"t\t"$3"\t"$4"\t"$5"\t"$6}' > ${TMPDIR}/${SAMPLE}.${i}keys
done

echo "--> Realigning mates of paired-reads with BLAT ..."

if [ -f ${TMPDIR}/${REFNAME}.11.ooc ];
then
   echo ""
else
   ${BLAT} ${REFGENOME} /dev/null /dev/null -tileSize=11 -makeOoc=${TMPDIR}/${REFNAME}.11.ooc -repMatch=1024
fi

echo ""

${BLAT} ${REFGENOME} ${TMPDIR}/${SAMPLE}.firstbothMates.TLSAmate.sorted.forBlat.seqs.fa  -ooc=${TMPDIR}/${REFNAME}.11.ooc ${TMPDIR}/${SAMPLE}.firstbothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out -out=blast8 
${BLAT} ${REFGENOME} ${TMPDIR}/${SAMPLE}.secondbothMates.TLSAmate.sorted.forBlat.seqs.fa  -ooc=${TMPDIR}/${REFNAME}.11.ooc ${TMPDIR}/${SAMPLE}.secondbothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out -out=blast8 

echo "--> Processing BLAT output ..."

for i in first second
do
	sort -k1,1 -k12,12nr ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out \
	| sort -u -k1,1 --merge  \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.bestHit

	cut -f 1 ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out \
	| sort \
	| uniq -c \
	| awk '{print $1"\t"$2}' \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.totalBlatHits

	join -1 1 -2 2 ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.bestHit ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.totalBlatHits \
	| join -1 1 -2 1 ${TMPDIR}/${SAMPLE}.${i}keys - \
	| python ${SP2TAB} \
	| cut -f2- \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combined

	perl ${REFORMATBLAT} ${TMPDIR}/${SAMPLE}.${i}Coord.regions ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combined ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregions
done

echo "--> Annotating reads and mates ..."

sort -k 2,2 ${L1INS} > ${TMPDIR}/${SAMPLE}.sorted

for i in first second
do
	awk '{print $0"\t"$1"\t"$2+500"\t"$3-500}' ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregions \
	| sort -k 23,23 \
	| join -1 23 -2 2 - ${TMPDIR}/${SAMPLE}.sorted \
	| python ${SP2TAB} \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea

	sort -k 1,1 ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat > ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.sorted
	sort -k 5,5 ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea \
	| join -1 5 -2 1 - ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.sorted \
	| python ${SP2TAB} \
	| sort -k 3,5 \
	| awk 'BEGIN { OFS="\t" } {print $3,$4,$5,$1,$28,$29,$30,$31,$32,$33,$34,$35,$36,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$27,$23,$2,$24}' \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted

done

echo "--> Filtering down the reads based on alignment length (>=${MINMATCHLEN}bp) and futher annotation ..."
for i in first second
do
	awk '$16>='$MINMATCHLEN ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted_AL

	cat ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted_AL \
	| awk '{ if ($19>1)             {lseq=substr($13,1,$19-1)}  else {lseq="0"}; \
			 if ($20<length($13)) {rseq=substr($13,$20)}      else {rseq="0"}; \
			 print $0 "\t" lseq "\t" rseq }' \
	| sort -k1,3 -k10,11 -k6,7 -k14,14 -k21,21 \
	> ${TMPDIR}/${SAMPLE}.${i}bothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted_AL_mismatch
done


cat ${TMPDIR}/${SAMPLE}.firstbothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted_AL_mismatch ${TMPDIR}/${SAMPLE}.secondbothMates.TLSAmate.sorted.forBlat.seqs.fa.blat.out.combinedWregionsWtea.reformatted_AL_mismatch \
| sort \
| uniq -c \
> ${TMPDIR}/${SAMPLE}.merged

echo "--> Preparing the final output ..."


cat ${TMPDIR}/${SAMPLE}.merged | awk 'BEGIN { OFS="\t" }{$1="";print $2"_"$3"_"$4,$0}' \
	| awk 'BEGIN { OFS="\t" } {if ($23>$22) {print $0} else {t=$23; $23=$22; $22=t; print $0}} '\
	| sort -k1,1 -k15,15 -k22,23n \
	| awk 'BEGIN{ OFS="\t"; id=""}\
		   { if ($1==id && $15==blatChr && $22<=max) \
		   		{ if ($23>max) max=$23; sum++; numMatchSum+=$26; if ($31$32~/AAAAAA|TTTTTT/) polyATcount++ } \
		   	else \
		   		{ if (sum>=4) {print id,blatChr,min,max,sum,numMatchSum/sum,polyATcount,$27 } id=$1; sum=1; min=$22; max=$23; blatChr=$15; numMatchSum=$26; if ($31$32~/AAAAAA|TTTTTT/) polyATcount=1; \
		   		else polyATcount=0 } } \
		   END{ if (sum>=4) { print id,blatChr,min,max,sum,numMatchSum/sum,polyATcount,$27 } }' \
	| awk 'BEGIN { OFS="\t" } {if (($6<=3)&&($5<30)) print $0}'\
	| sed 's/_/	/g'\
	| awk 'BEGIN { OFS="\t" } {if ($1!=$4) print $0}' \
	| awk '{print $1"\t"$2+500"\t"$3-500"\t"$0}' \
	| sort -n \
	> ${TMPDIR}/${SAMPLE}.passed_filters

awk 'BEGIN {print "TargetChr\tTargetStart\tTargetEnd\tSourceChr\tSourceStart\tSourceEnd\tSumOfReads\tAverageUniq\tReadsWith6A/Ts\tTSD\tSourceSize"}' > $OUTPUT
cat ${TMPDIR}/${SAMPLE}.passed_filters |awk 'BEGIN { OFS="\t" } {print $1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$9-$8}' >> $OUTPUT

NUMOFL1=$(wc -l $4 | awk '{print $1}')
NUMOFTRANS=$(wc -l $5 | awk '{print $1}')
NUMOFFINAL=$(wc -l ${TMPDIR}/${SAMPLE}.passed_filters | awk '{print $1}')

NUMOFL1=$(wc -l $4 | awk '{print $1}')
NUMOFTRANS=$(wc -l $5 | awk '{print $1}')
NUMOFFINAL=$(wc -l ${TMPDIR}/${SAMPLE}.passed_filters | awk '{print $1}')


echo "--> TIGER finished ..."

printf "\nFrom $NUMOFL1 L1 elements in total and $NUMOFTRANS translocations in total, $NUMOFFINAL L1-mediated sequence transductions were detected.\n"

#adding the hi-conf filtering
while true
do
	echo -n "Do you want to continue with the high-confidence filtering? (yes/no) "
	read input
	if [ "$input" == "yes" ]
	then
		echo -n "Enter the path to the corresponding Repeat Masker file in bed format [ENTER]: "
		read repmask
		awk 'BEGIN {print "TargetChr\tTargetStart\tTargetEnd\tSourceChr\tSourceStart\tSourceEnd\tSumOfReads\tAverageUniq\tReadsWith6A/Ts\tTSD\tSourceSize"}' > ${OUTPUT}_noRepMask
		${BEDTOOLS} intersect -a ${TMPDIR}/${SAMPLE}.passed_filters -b ${repmask} -v > ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask.bed
		cat ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask.bed >> ${OUTPUT}_noRepMask
		NUMOFRM=$(wc -l ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask.bed | awk '{print $1}')
		printf "\nDone filtering for repeats: $NUMOFRM transductions remained.\n"
		while true
		do
			echo -n "Do you have segmental duplication file in bed format available? (yes/no) "
			read seginput
			if [ "$seginput" == "yes" ]
			then
				echo -n "Enter the path of the corresponding segmental duplication file in bed format [ENTER]: "
				read segdups
				awk 'BEGIN {print "TargetChr\tTargetStart\tTargetEnd\tSourceChr\tSourceStart\tSourceEnd\tSumOfReads\tAverageUniq\tReadsWith6A/Ts\tTSD\tSourceSize"}' > ${OUTPUT}_noRepMask_noSegDups
				${BEDTOOLS} intersect -a ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask.bed -b ${segdups} -v > ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask_noSegDups.bed
				cat ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask_noSegDups.bed >> ${OUTPUT}_noRepMask_noSegDups
				NUMOFSG=$(wc -l ${TMPDIR}/${SAMPLE}.passed_filters_noRepMask_noSegDups.bed | awk '{print $1}')
				printf "\nDone filtering for segmental duplications: $NUMOFSG transductions remained.\n"
				echo "--> Finished with the high-confidence filtering ..."
				exit 1	

			elif [ "$seginput" == "no" ]
			then
				echo "No filtering based on presence of segmental duplications - quitting ..."
				echo "--> Finished with the high-confidence filtering ..."
				exit 1
			fi
		done	


	elif [ "$input" == "no" ]
	then
		echo "No filtering based on presence of reference repeats - quitting ..."
		exit 1
	fi
done
