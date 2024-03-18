#!/bin/bash

${PRIMERSEARCH} \
-seqall ${TRANSCRIPTOME} \
-infile <(cut -d$'\t' -f2,3,4 ${LVL}_targets.txt) \
-mismatchpercent ${PMISMATCH} \
-snucleotide1 \
-outfile ${LVL}_tsm.primersearch
