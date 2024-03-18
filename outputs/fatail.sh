#!/bin/bash

PTH=$(pwd)

# Add secondary structure info
# but automagically from the R execution
# + properly format the FASTA part
read -r NAME FOR_EXP REV_EXP <<< "$(echo -e "${GN}" | cut -d$'\t' -f1,2,3)"
for MET in $(ls "$PTH"/*"${TRANSCRIPT}".txt); do
	FOR_OBS=($(grep PRIMER_LEFT_._SEQUENCE $MET | cut -d= -f2))
	REV_OBS=($(grep PRIMER_RIGHT_._SEQUENCE $MET | cut -d= -f2))
	FOR_NUM=($(grep -n PRIMER_LEFT_._SEQUENCE $MET | cut -d: -f1))
	FOR_END=($(grep -n PAIR_._PRODUCT_TM $MET | cut -d: -f1))
	ALN="$PTH"/${NAME}_aln.txt
	for primer in ${!FOR_OBS[*]}; do
		[ "${FOR_OBS[primer]}" == "$FOR_EXP" ] && [ "${REV_OBS[primer]}" == "$REV_EXP" ] && sed -n "${FOR_NUM[primer]},${FOR_END[((primer))]}p" $MET | sed 's/\\n/\n/g' >> $ALN;
	done
	POOT=$(grep -c PRIMER3 $ALN)
	for pos in $(seq 0 $(($POOT-1))); do
		GRAPPENDIX=($(grep -n FASTA $ALN | cut -d: -f1))
		P3=($(grep -n PRIMER3 $ALN | cut -d: -f1))
		sed -i "${GRAPPENDIX[pos]}","${P3[pos]}"'s#\(^ENST.*\t\)#\n>\1\n#g' $ALN;
	done
done
