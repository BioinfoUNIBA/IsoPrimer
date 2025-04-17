#!/bin/bash
[ -e "${LVL}_psrch_resh.txt" ] && rm ${LVL}_psrch_resh.txt && touch ${LVL}_psrch_resh.txt

INPUT=${LVL}_tsm.primersearch
RESULT=${LVL}_psrch_resh.txt

IDX=($(grep -n Primer $INPUT | cut -d: -f1))
TEST=($(seq "${#IDX[*]}" | xargs))
IDX["${#IDX[*]}"]=$(($(wc -l $INPUT | cut -d' ' -f1)+1))

MAXLENGTH=$(($(grep RANGE rock.bakup | cut -d$'-' -f2)+25))

for e in ${TEST[*]}; do
	START=${IDX[(($e-1))]}
	STOP=$((${IDX[$e]}-1))
	MATCH=$(sed -n "${START},${STOP}p" $INPUT | grep 'Sequence' | sed 's/^.*: .*;\(.*\)$/\1/')
	TMATCH=($(sed -n "$(($START+1)),${STOP}p" $INPUT | grep 'Sequence' | sed 's/^.*: \(.*\);.*/\1/'))
	GENE=$(sed -n "${START},${STOP}p" $INPUT | grep Primer | cut -d' ' -f3)
	OUT=($(grep "length: [0-9]* bp" <(sed -n "${START},${STOP}p" $INPUT) | cut -d$' ' -f3))
	LONG=$(for k in ${!OUT[*]}; do j=${OUT[$k]} && (( $j > $MAXLENGTH )) && echo "${TMATCH[$k]}-$j"; done)
	ALL=$(for k in ${!OUT[*]}; do j=${OUT[$k]} && echo "${TMATCH[$k]}-$j"; done)
	echo -e "$GENE\t$(echo $MATCH)\t$(echo $LONG)\t$(echo $ALL)" >> $RESULT;
done
