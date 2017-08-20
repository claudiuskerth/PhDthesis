#!/bin/bash

EX=$1

for SFS in bootstrap/PAR/*dadi.corr;
do
	L=`tail -n 2 $SFS | head -n 1 | gawk '{for(i=1;i<=NF;i++)s+=$i;printf("%.0f\n", s)}'`; # get total sequence length from SFS
	head -n 1 $EX | gawk -v L=$L '{OFS="\t"}{$1="PAR";$3=L; print}' > $SFS.stair; # print custom header line
	echo -n "	" >> $SFS.stair; # this should print a TAB
	tail -n 2 $SFS | head -n 1 | sed -r 's/\s+/\t/g' | cut -f 2-19  >> $SFS.stair; # print SFS without 0th frequency class (monomorphic)
done
