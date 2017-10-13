#!/bin/bash

# simulate 200,000 reads of 50 bp length from the Heliconius reference;
# without sequencing error;
# with ploymorphism rate of 0.01, 20% of which are indels;
# indels are up to 20 bp long, length is randomly drawn from uniform distribution
# put read info in fastq header
#mason illumina -s 99 -N 200000 -o Mason_Hmel_sim.fq -sq -i -vv -n 50 -pi 0 -pd 0 -pmm 0 -pmmb 0 -pmme 0 -nN -hn 1 -hs 0.008 -hi 0.002 -hM 20 ../Hmel1-1_primary_contigs.fa 2> Mason_sim.log >> Mason_sim.log

# make fastq header wgsim_eval.pl compatible
#cat Mason_Hmel_sim.fq | perl -ne'unless (($.-1)%4==0){print; next;} /contig=(\w+).*orig_begin=(\d+).*orig_end=(\d+).*snps=(\d+).*indels=(\d+)/; print "\@$1_$2_$3_$4_$5\n";' > Mason_Hmel_sim_wgsimCompat.fq

###################################################################################
# Stampy
###################################################################################
# filter(){
# 	gawk '/^\@/ || ((and($2, 4)==0 || and($2, 8)==0) && $5>=4)' 
# 	# print either header line or if 4th or 8th bit set in SAM flag and mapping quality ge 4
# }
# export -f filter
# 
# bamsort(){
# 	samtools view -uhS - | samtools sort - Hme_sim_stampy
# }
# export -f bamsort
# 
# /usr/bin/time -o stampy.time stampy -g ../Hmel1-1_primary_contigs -h ../Hmel1-1_primary_contigs -t 11 --logfile Mason_Hme_sim_stampy.log -M Mason_Hmel_sim_wgsimCompat.fq | filter | tee Mason_Hme_sim_stampy.sam | bamsort
# 
# mutt -s "stampy mapping completed" claudius@huluvu < stampy.time 

# generate ROC curve
# mapping is correct if within 50 bp of true origin
# Note, that I needed to change the `wgsim_eval.pl` script slightly:
# On line 62, I needed to make the first pattern match for the chromosome/contig name non-greedy.
#../wgsim/wgsim_eval.pl alneval -ag 50 Mason_Hme_sim_stampy.sam > Mason_Hme_sim_stampy.roc

###################################################################################
# Razers3
###################################################################################
# 1
# this parameterisation should only report unique matches, but clearly it did not
# /usr/bin/time razers3 -fl swift -tc 10 -i 80 -rr 100 --unique -o razers3/Mason_sim_Hmel_razers3.sam -v ../Hmel1-1_primary_contigs.fa Mason_Hmel_sim_wgsimCompat.fq
# 
# ../wgsim/wgsim_eval.pl alneval -ag 50 razers3/Mason_sim_Hmel_razers3.sam > razers3/Mason_Hmel_sim_razers3.roc
# 
# mutt -s "ROC curve for razers3 completed" claudius@huluvu < nohup.out

# 2
# instead of '--unique', I just issue '-pa 1' that should purge reads that have more than one best match,
# i. e. two or more matches within the allowed edit distance (-i).
# /usr/bin/time razers3 -fl pigeonhole -tc 10 -i 80 -rr 100 -pa -o razers3/Mason_sim_Hmel_razers3_v2.sam -v ../Hmel1-1_primary_contigs.fa Mason_Hmel_sim_wgsimCompat.fq
# 
# ../wgsim/wgsim_eval.pl alneval -apg 50 razers3/Mason_sim_Hmel_razers3_v2.sam > razers3/Mason_Hmel_sim_razers3.roc 2> razers3/wrong.alingments_v2
# 
# mutt -s "ROC curve for razers3 completed" claudius@huluvu < nohup.out
# I needed to abort this run after 7h41m on 10 CPU's and a maximum memory usage of 91GiB as reported by `time`.
# It seems that an identity threshold of 80% is too low to allow razers3 to find all possible matches for the 200,000
# in the Heliconius reference sequence.

###################################################################################
# Bowtie2
###################################################################################

#bowtie2-build -f ../Hmel1-1_primary_contigs.fa Hmel1-1_primary_contigs_bowtie

# 1
# bowtie2 in end-to-end alignment mode
# /usr/bin/time bowtie2 -x ./bowtie2/Hmel1-1_primary_contigs_bowtie -U Mason_Hmel_sim_wgsimCompat.fq -S ./bowtie2/Mason_sim_Hmel_bowtie2.sam \
# 	-q --phred33 --very-sensitive --end-to-end -t -p 10
# 
# ../wgsim/wgsim_eval.pl alneval -apg 50 ./bowtie2/Mason_sim_Hmel_bowtie2.sam > ./bowtie2/Mason_sim_Hmel_bowtie2.roc 2> ./bowtie2/Mason_sim_Hmel_bowtie2.wrongAlign.sam
# 
# mutt -s "ROC curve for bowtie2 completed" claudius@huluvu < nohup.out

# 2
# bowtie2 in local alignment mode
# /usr/bin/time bowtie2 -x ./bowtie2/Hmel1-1_primary_contigs_bowtie -U Mason_Hmel_sim_wgsimCompat.fq -S ./bowtie2/Mason_sim_Hmel_bowtie2.sam \
# 	-q --phred33 --very-sensitive-local --local -t -p 10
# 
# ../wgsim/wgsim_eval.pl alneval -apg 50 ./bowtie2/Mason_sim_Hmel_bowtie2.sam > ./bowtie2/Mason_sim_Hmel_bowtie2.roc 2> ./bowtie2/Mason_sim_Hmel_bowtie2.wrongAlign.sam
# 
# mutt -s "ROC curve for bowtie2 completed" claudius@huluvu < nohup.out

###################################################################################
# Yara
###################################################################################

# /usr/bin/time yara_indexer -xp Hmel1-1_primary_contigs_yara -t ./yara ../Hmel1-1_primary_contigs.fa 
# 
# mutt -s "Yara genome index completed" claudius@huluvu < nohup.out

###################################################################################
# BWA
###################################################################################

# /usr/bin/time bwa index -p Hmel1-1_primary_contigs ../Hmel1-1_primary_contigs.fa
# 
# mutt -s "BWA genome index completed" claudius@huluvu < nohup.out

# 1
# default setting
# /usr/bin/time bwa mem -t 10 -T 1 ./bwa/Hmel1-1_primary_contigs_bwa Mason_Hmel_sim_wgsimCompat.fq > ./bwa/Mason_sim_Hmel_bwa_MEM.sam
# 
# mutt -s "BWA mapping completed" claudius@huluvu < nohup.out

# ../wgsim/wgsim_eval.pl alneval -apg 50 ./bwa/Mason_sim_Hmel_bwa_MEM.sam > ./bwa/Mason_sim_Hmel_bwa_MEM.roc 2> ./bwa/Mason_sim_Hmel_bwa.wrongAlign.sam


# 2
# more sensitive settings, i. e. with seed length 17, re-seeding if MEM>17*1.2 (=20.4) and mismatch penalty is 3,
# which together with a match score of 1 should have a target identity of 98.8% (or target divergence of 1.2%).
# I simulated the reads with a polymorphism frequency of 1%.
# /usr/bin/time bwa mem -t 10 -T 1 -k 17 -r 1.2 -B 3 ./bwa/Hmel1-1_primary_contigs_bwa Mason_Hmel_sim_wgsimCompat.fq > ./bwa/Mason_sim_Hmel_bwa_MEM.sam
# 
# mutt -s "BWA mapping completed" claudius@huluvu < nohup.out

#../wgsim/wgsim_eval.pl alneval -apg 50 ./bwa/Mason_sim_Hmel_bwa_MEM.sam > ./bwa/Mason_sim_Hmel_bwa_MEM.roc 2> ./bwa/Mason_sim_Hmel_bwa.wrongAlign.sam


###################################################################################
# SMALT
###################################################################################

# creating a hash index of the Heliconius genome reference sequence with kmer length 13 and 
# step/incremenet size 5
# /usr/bin/time smalt index -k 13 -s 5 ../Hmel1-1_primary_contigs ../Hmel1-1_primary_contigs.fa

# mapping with smalt on 10 cores, reads with multiple best mappings are reported as unmapped
# and default SW penalty scores
# /usr/bin/time smalt map -F fastq  -n 10 -r -1 -S match=1,subst=-2,gapopen=-4,gapext=-3 -o ./smalt/Mason_sim_Hmel_smalt.sam ../Hmel1-1_primary_contigs Mason_Hmel_sim_wgsimCompat.fq 
# 
# mutt -s "Smalt mapping completed" claudius@huluvu < nohup.out
# 
# ../wgsim/wgsim_eval.pl alneval -apg 50 ./smalt/Mason_sim_Hmel_smalt.sam > ./smalt/Mason_sim_Hmel_smalt.roc 2> ./smalt/Mason_sim_Hmel_smalt.wrongAlign.sam


###################################################################################
# Yara
###################################################################################

# Yara's genome indexer only allows up to 256 contigs. So I needed to concatenate the >11,000 contigs in the
# Heliconius genome reference with concat_ref.pl. 
# Since the coordinates and contig names are different in the concatenated reference,
# I need to redo the simulation of reads.

# simulate 200,000 reads of 50 bp length from the Heliconius reference;
# without sequencing error;
# with ploymorphism rate of 0.01, 20% of which are indels;
# indels are up to 20 bp long, length is randomly drawn from uniform distribution
# put read info in fastq header
# mason illumina -s 99 -N 200000 -o Mason_Hmel_sim_ref_concat.fq \
# 	-sq -i -vv -n 50 -pi 0 -pd 0 -pmm 0 -pmmb 0 -pmme 0 -nN -hn 1 -hs 0.008 -hi 0.002 -hM 20 ../Hmel1-1_primary_contigs_concat.fa 2> Mason_sim_ref_concat.log >> Mason_sim_ref_concat.log
# 
# # make fastq header wgsim_eval.pl compatible
# cat Mason_Hmel_sim_ref_concat.fq | \
#        	perl -ne'unless (($.-1)%4==0){print; next;} /contig=(\w+).*orig_begin=(\d+).*orig_end=(\d+).*snps=(\d+).*indels=(\d+)/; print "\@$1_$2_$3_$4_$5\n";' > Mason_Hmel_sim_ref_concat_wgsimCompat.fq
# 
# mutt -s "mason simulation finished" claudius@huluvu < nohup.out

# # create index
# /usr/bin/time yara_indexer ../Hmel1-1_primary_contigs_concat.fa

# run mapper with 10% maximal edit distance on 10 cores an 20,000 reads batch size
# /usr/bin/time yara_mapper -o ./yara/Mason_sim_Hmel_ref_concat_yara.sam -e 10 -t 10 -r 20000 ../Hmel1-1_primary_contigs_concat.fa Mason_Hmel_sim_ref_concat_wgsimCompat.fastq

# mutt -s "yara mapping completed" claudius@huluvu < nohup.out
# 
# ../wgsim/wgsim_eval.pl alneval -apg 50 ./yara/Mason_sim_Hmel_ref_concat_yara.sam > ./yara/Mason_sim_Hmel_ref_concat_yara.roc 2> ./yara/Mason_sim_Hmel_ref_concat_yara.wrongAlign.sam




