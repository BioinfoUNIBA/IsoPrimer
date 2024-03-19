#!/bin/bash

# Create a fresh output file
OUT=kalcounts.tsv
TMP=tmp.txt
[ -e "$OUT" ] && rm $OUT; touch $OUT
[ -e "$TMP" ] && rm $TMP

# Retrieve counts' paths and compile a table
KAL=$(ls $(pwd)/*/abundance.tsv)
# Luckily transcripts appear in the same order across samples
for e in $KAL; do
	paste $OUT <(cut -d$'\t' -f5 $e) > $TMP && mv $TMP $OUT;
done

# Add an header
DIR=$(dirname $KAL)
SAM=$(basename -a $DIR)
tail -n +2 $OUT > $TMP
echo $SAM | tr ' ' '\t' > $OUT && cat $TMP >> $OUT
sed -i "1,11s/\n/\t/g" $OUT

# Add row names (transcript IDs)
cat $OUT > $TMP
paste <(cut -d$'.' -f1 "${KAL%%$'\n'*}" | sed "1s/.*/t_name/") $TMP > $OUT
sed -i "s/\t\+/\t/g" $OUT
rm $TMP
