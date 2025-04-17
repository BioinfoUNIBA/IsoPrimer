#/bin/bash

# This file is the main launcher of IsoPrimer
# the options required should be inserted in the dedicated lines as follows:
# number of threads
# path to the Kallisto executable or 'custom' if the transcscript quantitation is manually provided
# path to the Primer3 executable
# path to the EMBOSS PrimerSearch executable
# Expression threshold percentage 
# 	(e.g. 50= consider isos more expressed than the median of the variants TPM distribution
# 	and 75= isoforms in the upper quartile are considered expressed)
# PrimerSearch mismatch percentage allowed
# path to the transcriptome (in .fa format)
# path to the genomic annotation
# path to the genome (in .fa format)
# debug flag (T/F)

nohup Rscript IP.R \
	4 \
	kallisto \
	primer3_core \
	primersearch \
	50 \
	20 \
	./test/test_transcriptome.fa \
	./test/test_annotation.gtf \
	./test/test_genome.fa \
	F &
