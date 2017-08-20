#!/bin/bash

EX=$1

for SFS in bootstrap/PAR/*folded;
do
	L=`gawk '{for(i=1;i<=NF;i++)s+=$i;printf("%.0f\n", s)}' $SFS`; # get total sequence length from SFS
	head -n 1 $EX | gawk -v L=$L '{OFS="\t"}{$1="PAR";$3=L; print}' > $SFS.stair; # print custom header line
	echo -n "	" >> $SFS.stair; # this should print a TAB
	cut -f 2- $SFS >> $SFS.stair; # print SFS without 0th frequency class (monomorphic)
done
