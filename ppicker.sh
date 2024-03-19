#!/bin/bash

# Make a temporary copy of the standard boulder.io file
cp -f rock.bakup ${LVL}_rock.temp
# Make a primer list, flush previous versions
[ -f "${LVL}"_primers.txt ] && rm ${LVL}_primers.txt && touch ${LVL}_primers.txt

# A loop will feed primer3 all transcript sequences ∀ gene from this path
REF=${LVL}_transcript_seqs_and_junctions.txt
GENE=$(head -n 1 $REF | cut -d$'\t' -f1)

grep $GENE $REF > ${LVL}_sub.temp 
sed -i "s/\t\+/\t/g" ${LVL}_sub.temp
while IFS= read -r TRANSCRIPT; do
	read -r TNAME SEQ SJS <<< "$(echo -e "$TRANSCRIPT" | cut -d$'\t' -f2,3,4)"
	sed -i "s/\(SEQUENCE_ID=\).*/\1$TNAME/" ${LVL}_rock.temp
	sed -i "s/\(SEQUENCE_TEMPLATE=\).*/\1$SEQ/" ${LVL}_rock.temp
	[ -z "$SJS" ] && SJS=0
	for J in $SJS; do
		[ "$J" == 0 ] && sed -i "/SEQUENCE_OVERLAP_JUNCTION_LIST=/d" ${LVL}_rock.temp || sed -i "s/\(SEQUENCE_OVERLAP_JUNCTION_LIST=\).*/\1$J/" ${LVL}_rock.temp
		${PRIMER3} < ${LVL}_rock.temp > outputs/${LVL}_${TNAME}.txt
		for e in $(seq 0 4); do
			PRIMERS=$(grep ${e}_SEQUENCE outputs/${LVL}_${TNAME}.txt | cut -d'=' -f2)
			TMS=$(grep ${e}_TM outputs/${LVL}_${TNAME}.txt | cut -d'=' -f2)
			PSZ=$(grep ${e}_PRODUCT_SIZE outputs/${LVL}_${TNAME}.txt | cut -d'=' -f2)
			[ -z "$PRIMERS" ] && PRIMERS=$(echo -e "NA\tNA") || PRIMERS=$(printf "%b\t%b" $PRIMERS)
			[ -z "$TMS" ] && TMS=$(echo -e "NA\tNA") || TMS=$(LC_NUMERIC="en_US.UTF-8" printf "%.0f\t%.0f" $TMS)
			[ -z "$PSZ" ] && PSZ=$(echo -e "NA")
			echo -e "$GENE\t${TNAME}_${e}\t$PRIMERS\t$TMS\t$PSZ" >> ${LVL}_primers.txt;
		done
		cat outputs/${LVL}_${TNAME}.txt >> outputs/${LVL}_${TNAME}.bakup
	done
	mv outputs/${LVL}_${TNAME}.bakup outputs/${LVL}_${TNAME}.txt
done < ${LVL}_sub.temp

# Cleanup
sed -i "s/ \+/\t/g" ${LVL}_primers.txt
sed -i "s/\t\+/\t/g" ${LVL}_primers.txt
