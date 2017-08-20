#!/bin/bash

#$ -P molecol
#$ -q molecol.q

# request max amount of memory to be used for each task
# -l mem=12G
# Job name, used to name output files and in the queue list
#$ -N PAR_realSFS
# Execute the job from the directory specified 
#$ -wd /fastdata/bop08ck
# Request 48 cores in an OpenMP environment 
#$ -pe openmp 48
# request maximum wall clock runtime for each task
#$ -l h_rt=07:59:59
# send an email at the start and end of each task, or when aborts
#$ -m abe
#$ -M c.kerth@sheffield.ac.uk

# tell the SGE where to find your executable
export PATH=/home/$USER/bin:$PATH

# Set the OPENMP_NUM_THREADS environment variable to 32
export OMP_NUM_THREADS=48

# # change into the directory where the input files are
# cd /fastdata/$USER

LOG=60_to_199_1DSFS_PAR.log

# run realSFS on bootstrap replicates of SAF files

for ID in `seq 178 199 | perl -ne'printf("%03d\n", $_)'`;
do
	start=`date +%s`;
	SAF=minInd9_overlapping/SAF/bootstrap/PAR/$ID*idx;
	SFS=`echo $SAF | sed 's/SAF/SFS/' | sed 's/saf\.idx/sfs/'`; 
	realSFS -P 48 -maxIter 50000 -tole 1e-6 -m 0 $SAF 2>/dev/null > $SFS;
	end=`date +%s`;
	diff=`echo "scale=4; ($end-$start)/60" | bc`;
	echo -n "Finished replicate $ID.  " >> $LOG;
	printf "It took %.4f minutes to finish.\n" $diff >> $LOG;
done
