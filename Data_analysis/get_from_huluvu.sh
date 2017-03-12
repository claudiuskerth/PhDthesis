#!/bin/bash

HOME="/Users/Claudius/Documents/PhD/THESIS/kks32/LaTeX/Data_analysis/"
cd $HOME

cd /Users/Claudius/Documents/PhD/THESIS/kks32/LaTeX/Data_analysis/SNP-indel-calling/dadi
rsync -avz -e'ssh -l claudius' huluvu.shef.ac.uk:/data3/claudius/Big_Data/DADI/dadiExercises .
#
# sync assembly.sh, the bash script that contains all commands and comments for the assembly of
# the standard RAD data set (Big Data), i. e. it combines de novo assembly and genotype likelihood
# estimation with ANGSD:
cd /Users/Claudius/Documents/PhD/THESIS/kks32/LaTeX/Data_analysis/SNP-indel-calling/
rsync -avz -e'ssh -l claudius' huluvu.shef.ac.uk:/data3/claudius/Big_Data/assembly.sh .
#
# get Rmarkdown notebook for PCA and all its dependencies
# I have create a rsync-filters file (in rsync jargon called 'merge' file or 'dir-merge' file) in
# /data3/claudius/Big_Data/ANGSD on huluvu, if you change it, run the rsync command with --list-only first
# to see what files will be synced and how large they are
cd /Users/Claudius/Documents/PhD/THESIS/kks32/LaTeX/Data_analysis/SNP-indel-calling
Claudius$ rsync -avzh --partial --list-only --stats --filter=": rsync-filters" --prune-empty-dirs --max-size=10M -e'ssh -l claudius' huluvu.shef.ac.uk:/data3/claudius/Big_Data/ANGSD .
