#!/bin/bash

rm \
	nohup.out \
	IP_Log.out \
	sj_seqs.tsv \
	mod_transcriptome.fa \
	validation_candidates.txt \
	genomic_amplimers.txt \
	genomic_mismatch.txt \
	val_cand_4genomcheck.txt \
	./*RData \
	./*primersearch \
	./[0-9]*_*txt \
	outputs/[0-9]*_*txt

rm \
	./quantification/KA_CountingOutput/kalcounts.tsv \
	./quantification/KA_index/index.idx

find ./quantification/KA_CountingOutput -name abundance.tsv -exec bash -c 'rm -rd $(dirname {})' \;
