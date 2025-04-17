#!/bin/bash

[ -e KA_index/index.idx ] || [ -e KA_CountingOutput/kalcounts.tsv ] || "$KALLISTO" index \
-i ./KA_index/index.idx \
--make-unique \
"$TRANSCRIPTOME"

FOLDERPATHS=$(while IFS= read -r SAMPLE; do realpath "$SAMPLE"; done < sample_list.txt)
[ -e KA_CountingOutput/kalcounts.tsv ] || for path in $FOLDERPATHS; do
	FOLDERNAME=$(basename "$path")
	KLST_OUTFOLDER=KA_CountingOutput/"$FOLDERNAME"
	[ -d "$KLST_OUTFOLDER" ] || mkdir "$KLST_OUTFOLDER"
	FASTAFILE_FOR=$(ls "$path"/*1.fastq.gz)
	FASTAFILE_REV=$(ls "$path"/*2.fastq.gz)
	"$KALLISTO" quant \
	-i KA_index/index.idx \
	--bias \
	--rf-stranded \
	-b 150 \
	-t "$THREADS" \
	"$FASTAFILE_FOR" \
	"$FASTAFILE_REV" \
	-o "$KLST_OUTFOLDER";
done

################
#              #
#  single_end  #
#              #
################
#
# The default settings set above launch kallisto
# in paired-end mode with the --rf-stranded option
# substitute lines 15-23 with the following uncommented block to launch
# kallisto in single-end mode.
#	FASTAFILE_FOR=$(ls $path/*fastq.gz)
#	"$KALLISTO" quant \
#	-i KA_index/index.idx \
#	--bias \
#	--single \
#	-l 100 \
#	-s 1 \
#	--single \
#	-b 150 \
#	-t "$THREADS" \
#	"$FASTAFILE_FOR" \
#	-o "$KLST_OUTFOLDER";
# Modify the parameters of the -l and -s options if necessary.
# See https://pachterlab.github.io/kallisto/manual for more information

########################
#                      #
#  implementing_nohup  #
#                      #
########################
#
# The for cycle above quantifies the transcripts of a sample at a time
# IF ENOUGH RESOURCES ARE AVAILABLE speed up the quantification by 
# increasing the number of threads and substituting lines 15-23 with:
#	nohup "$KALLISTO" quant \
#	-i KA_index/index.idx \
#	--bias \
#	--rf-stranded \
#	-b 150 \
#	-t "$THREADS" \
#	"$FASTAFILE_FOR" \
#	"$FASTAFILE_REV" \
#	-o "$KLST_OUTFOLDER" &
