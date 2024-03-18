#!/bin/bash

grep "^>" "${TRANSCRIPTOME}" > og_headers.txt

grep -q 'biotype' og_headers.txt &&
	paste -d_ \
		<(sed 's#>\(ENS[^ ]*\).*\(ENS[^ ]*\).*#\1#' og_headers.txt) \
		<(sed 's#>\(ENS[^ ]*\).*\(ENS[^ ]*\).*#\2#' og_headers.txt) \
		<(sed -e '/.*gene_symbol.*/!s/.*//;s#.*gene_symbol:\([^ $]*\).*#\1#' og_headers.txt) \
		<(sed -e '/.*gene_biotype.*/!s/.*//;s#.*gene_biotype:\([^ $]*\).*#\1#' og_headers.txt | tr '_' '-') \
		> mod_headers.txt ||
sed 's#[>|]\([^|]*\)|\([^|]*\)|\([^|]*\)|\([^|]*\)|\([^|]*\)|\([^|]*\)|\([^|]*\)|\([^|]*\)|.*#\1_\2_\5_\8#g' og_headers.txt > mod_headers.txt

Rscript -e "suppressPackageStartupMessages(library(Biostrings))
 tm <- Sys.getenv('TRANSCRIPTOME')
 seqs <- readDNAStringSet(tm)
 new_heads <- scan('mod_headers.txt', what='character', quiet=T)
 names(seqs) <- new_heads
 writeXStringSet(seqs, 'mod_transcriptome.fa')"

rm og_headers.txt mod_headers.txt
