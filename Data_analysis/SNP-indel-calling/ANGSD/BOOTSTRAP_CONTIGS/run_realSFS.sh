#!/bin/bash

# run realSFS on bootstrap replicates of SAF files

for ID in `seq 0 99 | perl -ne 'printf("%03d\n", $_)'`; 
do 
	start=`date +%s`
	# sleep 1
 	SAF=including_non-overlapping/SAF/bootstrap/PAR/$ID*idx;
 	SFS=`echo $SAF | sed 's/SAF/SFS/' | sed 's/saf\.idx/sfs/'`; 
#  # 	echo -n $SAF "   ";
#  # 	echo $SFS;
 	realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 $SAF 2>/dev/null > $SFS
	end=`date +%s`
	diff=`echo "scale=4; ($end-$start)/60" | bc`
	echo -n "Finished replicate $ID.  " 1>&2
	printf "It took %.4f minutes to finish.\n" $diff
done
