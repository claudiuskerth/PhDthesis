#!/bin/bash

# imulation of 200,000 50bp SE reads from the Heliconius melpomene reference sequence
# substitution rate 0.01, 20% of substitutions are indels with length randomly drawn 
# from the geometric distribution 0.9^(length-1)*0.1
# no sequencing error
# see ttps://github.com/lh3/wgsim
#./wgsim/wgsim -N200000 -e0 -150 -250 -S13 -r0.01 -R0.2 -X0.9 -d0 Hmel1-1_primary_contigs.fa Hme_sim_SE.fq /dev/null
#
##./wgsim/wgsim -N100000 -e0 -150 -250 -d200 -s50 -S17 -r0.01 -R0.2 -X0.9 Hmel1-1_primary_contigs.fa Hme_sim_PE.fq_1 Hme_sim_PE.fq_2 
#
##stampy -G Hmel1-1_primary_contigs Hmel1-1_primary_contigs.fa
##stampy -g Hmel1-1_primary_contigs -H Hmel1-1_primary_contigs
#
#filter(){
#	gawk '/^\@/ || ((and($2, 4)==0 || and($2, 8)==0) && $5>=4)' 
#	# print either header line or if 4th or 8th bit set in SAM flag and mapping quality ge 4
#}
#export -f filter
#
#bamsort(){
#	samtools view -uhS - | samtools sort - Hme_sim_SE_stampy
#}
#export -f bamsort
#
#/usr/bin/time -o stampy_SE.time stampy -g Hmel1-1_primary_contigs -h Hmel1-1_primary_contigs -t 11 --logfile Hme_sim_SE_stampy.log -M Hme_sim_SE.fq | filter | tee Hme_sim_SE_stampy.sam | bamsort
#
#mutt -s "stampy mapping completed" claudius@huluvu < stampy_SE.time 

./wgsim/wgsim_eval.pl alneval -g 20 

# The next line extracts the median length contigs from the Heliconius reference sequence. 
# I got the median length from assemblathon_stats.pl on 'Hmel1-1_primary_contigs.fa'.
# perl -e'$_ = <>; $h = $_; while(<>){ while(!/^>/){chomp; $seq .= $_; exit unless $_ = <>}; if( length($seq) == 11515 ){ print $h, "$seq\n"}; $h = $_; $seq = "";}' Hmel1-1_primary_contigs.fa > Hmel1-1_primary_contigs_red_ref.fa
