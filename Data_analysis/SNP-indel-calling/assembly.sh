#!/bin/bash
# 
# # # -----------------------------
# # # collapse all PCR duplicates
# # # -----------------------------
# # for se in data/*fq_1.gz;
# # do
# # 	pe=$(echo $se | sed 's/fq_1/fq_2/')
# # 	starcode -d 2 -t 12 -c -1 <(zcat $se) -2 <(zcat $pe) | gzip > nonredundant/`basename $pe .fq_2.gz`.collapsed.gz
# # done
# # # the -c switch should select one canonical sequence from each set of connected components in the graph, i. e. no clustering of graph, just raw networks
# # # this is equivalent to transitive clustering
# # # see Issue #13 on starcodes github repo
# # 
# # # -----------------------------
# # # count collapsed sequences
# # # -----------------------------
# # sum=0; 
# # for file in nonredundant/*collapsed.gz;
# # do 
# # 	echo -n "$file   " >> nonredundant/uniqseq.count;
# #        	count=`zcat $file | wc -l`;
# #        	sum=$[$sum + $count];
# #        	echo $count >> nonredundant/uniqseq.count;
# # done 
# # echo "sum   $sum" >> nonredundant/uniqseq.count
# # 
# # # ---------------------------------------------------------------------------------------
# # # combine collapsed single-end unique sequences from all individuals into one fasta file,
# # # including the individual ID in the fasta header
# # # ---------------------------------------------------------------------------------------
# # perl -e\
# # '@ARGV = map {"zcat $_ |"} @ARGV; while(<>){$n = $ARGV; ($n) = $n =~ m/redundant\/(.*)\.coll/; s/-+.*//; print ">Seq_", $., ";$n", "\n", $_;}' \
# # nonredundant/*collapsed.gz | \
# # gzip > nonredundant/totaluniqseq.SE.fasta.gz
# # # the perl command needs to modify @ARGV to read in gzipped files; $ARGV is automatically set to the name of the input file currently read in
# # # as perl iterates over @ARGV; $n is then set to the individual ID; SE and PE sequence are separated by dashes in the output from starcode, I only keep
# # # the SE sequences; I then give each sequence a unique number (line number from input lines), followed by the individual ID, followed by the sequence
# # 
# # perl -e\
# # '@ARGV = map {"zcat $_ |"} @ARGV; while(<>){$n = $ARGV; ($n) = $n =~ m/redundant\/(.*)\.coll/; s/-+.*//; print ">Seq_", $., ";$n", "\n", $_;}' \
# # nonredundant/ery*collapsed.gz | \
# # gzip > nonredundant/ery_totaluniqseq.SE.fasta.gz
# # #
# # perl -e\
# # '@ARGV = map {"zcat $_ |"} @ARGV; while(<>){$n = $ARGV; ($n) = $n =~ m/redundant\/(.*)\.coll/; s/-+.*//; print ">Seq_", $., ";$n", "\n", $_;}' \
# # nonredundant/par*collapsed.gz | \
# # gzip > nonredundant/par_totaluniqseq.SE.fasta.gz
# # 
# # # ------------------------------------------------
# # # cluster 'totaluniqseq.SE.fasta.gz' with vsearch
# # # ------------------------------------------------
# # cd VSEARCH
# # ln -s ../nonredundant/totaluniqseq.SE.fasta.gz .
# # nohup /usr/bin/time -vo time.log vsearch \
# # 	--cluster_fast totaluniqseq.SE.fasta.gz \
# # 	--id 0.8 --qmask dust \
# # 	--uc >(gzip > uclust.gz) \
# # 	--log totaluniqseq.SE.log &
# # cd .. 
# # # ----------------------------------------------------------------
# # # turn vsearch cluster output into 'rainbow cluster' output format
# # # ----------------------------------------------------------------
# # cd VSEARCH
# # zcat uclust.gz | mawk '!/^C/' | \
# # 	perl -ne'if(/^S/){@S=split; $s{$S[8]}=$S[1]; print;}elsif(/^H/){@H=split; $H[1]=$s{$H[9]}; print join("\t", @H), "\n"}' | \
# # 	cut -f 2,9 | sed -r 's/Seq_//' | sed 's/;.*//' | sort -gk2,2 | awk '{print $2 "\t" $1}' | \
# # 	paste - <(zcat ../nonredundant/*collapsed.gz | cut -f 1 | sed -r 's/-+/\t/') | \
# # 	sort -gk2,2 | \
# # 	mawk '$2=$2+1' | tr ' ' '\t' | \
# # 	gzip > rcluster.gz
# # # takes centroids and hits from uclust format output of vsearch (lines starting with a C contain cluster info); 
# # # perl line turns ordinal number of matching centroid for hit sequences (H) into cluster ID, see https://github.com/torognes/vsearch/issues/207
# # # takes cluster id and query id columns; strips everything away from query id except seq id number; sorts by seq id number; swaps cluster id and seq id around; 
# # # pastes in the unique sequence pairs from starcode output to create a TAB separated file with four columns: SEQ_ID<TAB>CLUSTER_ID<TAB>SE_SEQ<TAB>PE_SEQ 
# # # sorts by cluster ID
# # # adds 1 to 0-based cluster ID's from Vsearch; replaces space with TAB
# # cd ..
# # 
# # # -------------
# # # rainbow div
# # # -------------
# # cd RAINBOW
# # ln -s ../VSEARCH/rcluster.gz .
# # nohup /usr/bin/time -vo rainbow_div.time \
# # 	rainbow div -i <(gzip -dc rcluster.gz) -o >(gzip > rbdiv.out.gz) -k 2 -f 0.05 -K 10
# # cd ..
# # # if the minor sequence in a cluster has count of >= 2, then rainbow should split the cluster in two;
# # # if the minor sequence in a cluster has frequency >= 0.05, then rainbow should split the cluster; 
# # # both criteria are probably applied recursively to all split clusters until one of them is not met 
# # 
# # # ----------------------------------
# # # rainbow merge and contig assembly
# # # ----------------------------------
# # cd RAINBOW
# # /usr/bin/time -vo rb_merge.time rainbow merge -i <(zcat rbdiv.out.gz) -o >(gzip > rb.merge.asm.out.gz) -a -r 2
# # # require only two reads for contig assembly
# # cd ..
# #
# # # -----------------------------
# # # extract best contigs + read1
# # # -----------------------------
# # cd RAINBOW
# # zcat rb.merge.asm.out.gz | cat - <(echo "E") | sed 's/[0-9]*:[0-9]*://g' | mawk ' {
# # if (NR == 1) e=$2;
# # else if ($1 ~/E/ && lenp > len1) {c=c+1; print ">Contig_" e "\n" seq2 "NNNNNNNNNN" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
# # else if ($1 ~/E/ && lenp <= len1) {c=c+1; print ">Contig_" e "\n" seq1; seq1=0; seq2=0;lenp=0;e=$2;fclus=0;len1=0;freqp=0;lenf=0}
# # else if ($1 ~/C/) clus=$2;
# # else if ($1 ~/L/) len=$2;
# # else if ($1 ~/S/) seq=$2;
# # else if ($1 ~/N/) freq=$2;
# # else if ($1 ~/R/ && $0 ~/0/ && $0 !~/1/ && len > lenf) {seq1 = seq; fclus=clus;lenf=len}
# # else if ($1 ~/R/ && $0 ~/0/ && $0 ~/1/) {seq1 = seq; fclus=clus; len1=len}
# # else if ($1 ~/R/ && $0 ~!/0/ && freq > freqp && len >= lenp || $1 ~/R/ && $0 ~!/0/ && freq == freqp && len > lenp) {seq2 = seq; lenp = len; freqp=freq}
# # }' | \
# # perl -ne'if($. % 2){$h=$_}else{if(/CCTGCA$/){print $h, $_}}' \
# # > >(gzip > rainbow.fasta.gz)
# # cd ..
# # # from dDocent; In the output of `rainbow merge -a` every cluster starts with a line beginning with E and the cluster id. For each cluster several contigs (C) can
# # # be reported. Lines starting with N contain the number of sequences that went into the assembly of the contig (freq). Lines starting with an R contain sequence id's
# # # in the format <seq_id>:<0-based position of the sequence in the contig>:<[01]>. If the id string of a sequence ends with 1 after the second colon, then this is a paired-end
# # # sequence, if it ends with 0, then this is a single-end sequence and the seq_id before the first colon is identical to the cluster id. Usually there is a short
# # # single-end contig (RAD tag) and a longer paired-end contig. The awk code above prints this out in fasta format, the PE contig before the SE RAD tag, both separated
# # # by a string of N's. Sometimes PE reads overlap with the SE RAD tag and `rainbow merge -a` creates a contig including the SE RAD tag. If this is shorter than the
# # # PE contig, then this is printed out after the string of N's. Otherwise it is printed out on its own, ignoring the PE contig. Generally the longest PE contig is printed out.
# # # If two PE contigs from a cluster happen to be the same length, then the one with higher coverage is printed out.
# # # The perl line that follows makes sure that every sequence ends with "CCTGCA", the remainder of the SbfI restriction site. All SE reads that went into the analysis
# # # were filtered (or corrected) to start with this sequence. If it is missing at the end, then the `rainbow merge -a` assembly must be wrong.
# # 
# # 
# # # --------------------------------------------
# # # cluster rainbow assembly again with vsearch
# # # --------------------------------------------
# # cd VSEARCH
# # ln -s ../RAINBOW/rainbow.fasta.gz .
# # nohup /usr/bin/time -vo time.log vsearch \
# # 	--cluster_fast <(seqtk seq -r rainbow.fasta.gz) \
# # 	--id 0.90 --qmask dust \
# # 	--clusterout_id \
# # 	--clusterout_sort \
# # 	--centroids >(gzip > rainbow.vsearch.centroids.gz) \
# # 	--consout >(gzip > rainbow.vsearch.consensus.gz) \
# # 	--uc >(gzip > rainbow.vsearch.uclust.gz) \
# # 	--sizeout \
# # 	--log rainbow.fasta.vsearch.log &
# # cd .. 
# 
# 
# # # -------------
# # # READ MAPPING
# # # -------------
# # mkdir BOWTIE2 
# # cd BOWTIE2
# # ln -s ../VSEARCH/rainbow.vsearch.consensus.gz .
# # zcat rainbow.vsearch.consensus.gz | seqtk seq -l 1000 | perl -pe'if(/^>/){s/.*(Contig_\d+);.*/>$1/}else{s/^(.*)$/N$1N/}' | bgzip > Big_Data_ref.fa.gz
# # samtools faidx Big_Data_ref.fa.gz
# # cd ..
# # # the assembled RADome in fasta format has sequences broken over several lines, with 'seqtk' I put each reference sequence on a single line (no reference
# # # sequence is longer than 1000 bases); the perl command shortens the fasta headers only keeping the contig ID and adds an N at the beginning and end of 
# # # each sequence (according to Jon Puritz from dDOcent, Freebayes has a bug that can be circumvented by adding N's, see https://github.com/ekg/freebayes/issues/174); 
# # # finally the reference is bgzipped
# # 
# # 
# # # ---------------
# # # bowtie2 index
# # # ---------------
# # cd BOWTIE2
# # # need to create uncompressed reference fasta file for 'bowtie2-build'; it won't accept process substitution:
# # bgzip -dc Big_Data_ref.fa.gz > Big_Data_ref.fa
# # bowtie2-build  Big_Data_ref.fa  Big_Data_ref -f -o 4 -t 12 --threads 23
# # rm Big_Data_ref.fa
# # cd ..
# # # since huluvu has a lot of memory, I am slightly increasing the offrate and lookup table of the index, both should make searching the index faster
# # 
# # ----------------
# # bowtie2 mapping
# # ----------------
# # cd BOWTIE2
# # mkdir BAM
# # for SE in ../data/*fq_1.gz; do
# # 	PE=`echo $SE | sed 's/fq_1/fq_2/'`
# # 	ID=`basename $SE .fq_1.gz` 
# # 	bowtie2 -x Big_Data_ref \
# # 		-1 <(zcat $SE) -2 <(zcat $PE) \
# # 		-q \
# # 		--phred64 \
# # 		--very-sensitive \
# # 		--dpad 10 \
# # 		--gbar 4 \
# # 		--norc \
# # 		--end-to-end \
# # 		--np 10 \
# # 		-I 60 \
# # 		-X 800 \
# # 		--no-contain \
# # 		-t \
# # 		--no-unal \
# # 		--rg-id $ID \
# # 		--rg "SM:$ID" \
# # 		--omit-sec-seq \
# # 		-p 22 \
# # 		| samtools view -bq 1 - \
# # 		> BAM/$ID.bam;
# # done
# # cd ..
# # # running bowtie2 in end-to-end mode (i. e. no clipping of query sequences), "very sensitive" specifies a "seed" (kmer) length of 20
# # # that needs to match exactly, seeds are sampled from the query according to the function f(x) = 1 + 0.5 * sqrt(x), x is the read length
# # # with "--dpad 10" I am allowing gaps to be up to 10 bp long; "--gbar 4" disallows gaps within 4 bp of either end of the read; with
# # # "--norc" bowtie2 will not attempt to align unpaired reads against the revcomp reference strand - I generally don't expect any good alignments
# # # against the revcomp of the RADome, so this flag somewhat reduces unnecessary searches; "--np 10" sets the penalty for ambiguous 
# # # characters in an alignment - in the RADome SE RAD tags and PE contigs are separated by up to 10 N's, I don't want alignments across this gap;
# # # I specify a minimum fragment length of 60 (i. e. allowing some overlap between SE and PE reads) and a maximum fragment length of 800 (longer
# # # than the longest sequence in the RADome; I am supressing the output of SAM records for reads that failed to align anywhere; if secondary sequences
# # # are printed (which they shouldn't with this configuration), then omit the seq and qual string; "-p 22" runs bowtie2 on 22 cores; samtools view
# # # then filters out all SAM records with 0 mapping quality and writes out a compressed BAM file
# 
# # # -----------------------------------------------
# # # analyse distribution of SAM flags in BAM file
# # # ----------------------------------------------
# # samtools view par_34-9.bam | cut -f 2 | 
# # perl -ne'chomp; $flag{$_}++;END{while( ($k,$v)=each %flag){printf "%3d\t%6d\t%s", $k, $v, `samtools flags $k | cut -f3`;}}' | sort -gr -k2,2
# 
# # # -------------------------------------------------
# # # sort BAM files according to position in parallel
# # # -------------------------------------------------
# # parallel 'samtools sort -o {.}.sorted.bam {}' ::: ery*bam
# # parallel 'samtools sort -o {.}.sorted.bam {}' ::: par*bam
# # 
# # # -------------------------------------
# # # change RG tag: add SM, LB and PL tag
# # # -------------------------------------
# # cd BOWTIE2/BAM
# # rename -v 's/bam/bamx/' *sorted.bam
# # for f in *bamx; do 
# # 	n=`echo $f | sed 's/x$//'`;
# # 	samtools view -h $f | \
# # 	sed -r 's/^@RG\tID:(.*?)\t.*/@RG\tID:\1\tSM:\1\tLB:standard_RAD\tPL:GAIIx/' | \
# #        	samtools view -b > $n; 
# # 	rm -f $f;
# # done 
# # cd ../..
# # # 
# # # --------
# # # FILTERS
# # # --------
# # # -----------------------
# # # mismappings of SE reads
# # # -----------------------
# # cd BOWTIE2/BAM
# # samtools cat -o >(samtools view -f64 - | gawk '$4>2' | cut -f 3 | sort | uniq > exclude_contigs.SEposgt2) *bam
# # # With this command line I am creating a list of contig ID's where a SE read from any ind mapped to a position greater than 2 
# # # (remember, I padded the contigs with 1 N at beginning and end). A mapping position > 2 is a mismapping and indicates an
# # # incorrect assembly of the contig. The whole contig should therefore be taken out of any further analysis. 
# # # With the following command I create a list of all contigs that get at least one read mapped by any of the 36 individuals:
# # samtools cat *sorted.bam | samtools view | cut -f 3 | uniq | sort | uniq > Big_Data_Contigs_with_mappings.list
# # # The list contains 582,596 contigs. tHE WHOLE bIG dATA REFERENCE CONTAINS THIS MANY CONTIGS:
# # zcat ../Big_Data_ref.fa.gz | grep -c "^>"
# # THIS PRINTS 583,312. So almost all reference contigs got a read mapped.
# # # There are 4575 contigs with SE reads mapping beyond position 2. In order to create a positive list of contigs to include for ANGSD
# # # analysis, I have to subtract them from the list of all contigs that got reads mapped:
# # perl -ne'$H{$_}=1; END{open($IN, "Big_Data_Contigs_with_mappings.list"); while(<$IN>){print $_ unless exists $H{$_}}}' exclude_contigs.SEposgt2
# # # This perl command first reads in all contig ID's to exclude, then opens the list with all mapped contigs
# # # and prints each contig ID out if it is not contained in the hash of contig ID's to exclude. I thus get a positive list of
# # # contig ID's without SE mismappings.
# 
# 
# # # ---------------------------------
# # # low complexity sequence intervals
# # # ---------------------------------
# # # cd ..
# # bgzip -dc Big_Data_ref.fa.gz | dustmasker -out Big_Data_ref.dust_intervals -window 64 -level 20 -parse_seqids -outfmt acclist
# # # creates a list of intervals that where masked because of low complexity by dustmasker
# # # the intervals positions are 0-indexed and include the first and the last position
# # # BED format also uses 0-indexing of chr positions, but chrEND is the pos of the nt right after the feature (same as indexing in Python)
# # # In order to transform dust intervals to BED format, I therefore have to add 1 to the third column, i. e. chrEnd+1
# # # The following command creates a position sorted BED file from the dust intervals:
# # sed -r 's/>lcl\|//' Big_Data_ref.dust_intervals | \
# # 	sed 's/_/\t/' | \
# # 	sort -gk2,2 -k3,3 | \
# # 	perl -pe's/^(.*?)\t(\d)/$1_$2/' | \
# # 	gawk 'BEGIN{OFS="\t"}{$3+=1; print}' \
# # 	> Big_Data_ref.dust_intervals.sorted.bed
# # # I want to create the complement of these intervals with 'bedtools complement', i. e. create a positive list of intervals to include in the ANGSD analysis
# # # the following command creates a genome index file that is required for 'bedtools complement':
# # samtools view -h BAM/par_34-9.sorted.bam | \
# # 	grep "^@SQ" | \
# # 	sed 's/^.*SN://' | \
# # 	sed 's/LN://' | \
# # 	sed 's/_/\t/' | \
# # 	sort -gk2,2 | \
# # 	perl -pe's/^(.*?)\t(\d)/$1_$2/' \
# # 	> Big_Data_ref.idx
# # # The command takes the info from @SQ lines in the header of a BAM file that contains contig name as SN tag and contig length as LN tag.
# # # The contigs are then sorted by numerical id, similar to the position sorting of dust intervals during transformation to BED format.
# # # The following command creates a positive list of high-complexity intervals to include in ANGSD analysis
# # bedtools complement -i Big_Data_ref.dust_intervals.sorted.bed -g Big_Data_ref.idx > Big_Data_ref.dust_intervals.sorted.compl.bed
# # # This still needs to be transformed to ANGSD format:
# # gawk 'BEGIN{OFS="\t"}{$2+=1; print}' Big_Data_ref.dust_intervals.sorted.compl.bed > Big_Data_ref.dust_intervals.sorted.compl.angsd_sites
# # # ANGSD uses 1-indexing of positions and includes regStart and regEnd in the interval
# 
# 
# 
# # # -----------------------
# # # excessive depth filter
# # # -----------------------
# # # I would like to know the depth distribution over SE RAD tags. The coverage across PE contigs will necessarily be more dependent 
# # # on the length of the assembled contig.
# # cd /data3/claudius/Big_Data/BOWTIE2/BAM
# # getDepth(){
# # 	gawk '$4==2' | cut -f 3 | perl -ne'chomp; $H{$_}++;END{for $contig (keys %H){ print $contig, "\t", $H{$contig}, "\n";}'
# # }
# # parallel -j 10 "samtools view -f 64 {} | getDepth > {.}.SE_depths" ::: *sorted.bam &
# # # I need to put part of the code into a shell function in order to hide it from 'parallel' parsing. I am taking each individual BAM
# # # file, extract all SE reads from it, make sure each maps to the correct position (2). Then I take the contig ID's and count them in
# # # Perl hash. Finally I print out the coverage counts from each contig. This is done in parallel for each individual BAM file. The
# # # resulting files can be read into R for plotting, which I have done in Big_Data.Rmd.
# # # I decided that 25x coverage (~Q95) is a reasonable maxDepth threshold. I therefore want to exclude any contig from ANGSD analysis
# # # that has >25x coverage in any individual:
# # cat *depths | gawk '$2>25' | cut -f 1 | \
# # 	perl -ne'$H{$_}=1;END{open($IN, "Big_Data_Contigs.noSEgt2"); while(<$IN>){print $_ unless exists $H{$_};}}' | \
# # 	gawk '{print $1":"}' \
# # 	> Big_Data_Contigs.noSEgt2.depthLT26.angsd_regions
# # # This command takes all 36 depths files and filters out those lines with coverage greater 25. It then keeps only the contig ID's
# # # and reads them into a hash, thereby uniq'ing them, i. e. keeping only one rep from each contig ID. I then read in the contig ID's
# # # that passed the mismappings filter and print them out if they are not among the contigs with coverage >25x.
# #
# # # an alternative way to measure coverage at each position in a contig that got any reads mapped to.
# # bamtools coverage -in par_34-14.sorted.bam | gzip > par_34-14.sorted.cov.gz
# # bamtools coverage -in ery_30-15.sorted.bam | gzip > ery_30-15.sorted.cov.gz
# # # I am taking two arbitrarily chosen individuals from each population and determine the coverage for each position
# # # in R, I am then determining the coverage distribution and the 99th percentile of the coverage distribution for each individual
# # # this shall be the depth thresholds for each population above which a contig gets excluded from further analysis
# # # For ery I picked the 95th percentile, 19, as the threshold above which contigs get excluded, for par the 95th percentile is 18.
# # # see Big_Data.Rmd
# # zcat ery_30-15.sorted.cov.gz | gawk '$3>19' | cut -f 1 | uniq > ery_30-15.covGt19.exclude_contigs
# # zcat par_34-14.sorted.cov.gz | gawk '$3>18' | cut -f 1 | uniq > par_34-14.covGt18.exclude_contigs
# # Note, I opted for the first excessive depth filter based only on depth of SE RAD tags, see above.
# # 
# # #
# # # print the total length of the Big Data reference
# # #
# # cd ..
# # seqtk seq -l 1000 ../Big_Data_ref.fa.gz | gawk '{if(NR%2==0) sum+=length} END{print sum}'
# # cd BAM
# # # the total length is 98,732,423
# # 
# # zcat ery_30-15.sorted.cov.gz | wc -l
# # # prints out 23,261,159
# # # this should be the length of those contigs to which ery_30-15 had reads mapped to 
# # 
# # zcat ery_30-15.sorted.cov.gz | cut -f 1 | uniq | wc -l
# # # this prints 141,042
# # # so this many contigs from the Big Data reference got reads mapped to
# # # this can be verified by counting the number of contigs in the BAM file of this individual
# # samtools view ery_30-15.sorted.bam | cut -f 3 | uniq | wc -l
# 
# # # get length of longest contig in the Big Data reference
# # gawk 'NR%2==0' Big_Data_ref.fa | gawk '{print length}' | gawk 'BEGIN{max=0}{if($1>max)max=$1}END{print max}'
# # # 787
# 
# 
# # # ------------------
# # # INDEL realignment
# # # ------------------
# # Picard and GATK need java version 1.8. Select with $ update-alternatives --config
# # mkdir IndelRealignment
# # cd IndelRealignment
# # # need uncompressed reference sequence for GATK, a .fai (samtools faidx) and .dict file thatn can be created with PICARD:
# # java -jar $picard CreateSequenceDictionary R=Big_Data_ref.fa O=Big_Data_ref.dict
# # cd ..
# # ls *sorted.bam | parallel "samtools index {}"
# # cd BOWTIE2/BAM/IndelRealignment
# # #
# # # !!! This takes forever, Do not execute !!!
# # #
# # /usr/bin/time -vo RealignerTargetCreator.time java -Xmx60g -jar $gatk \
# # 	-T RealignerTargetCreator \
# # 	-R Big_Data_ref.fa \
# # 	-o forIndelRealigner.intervals \
# # 	-nt 12 \
# # 	--maxIntervalSize 800 \
# # 	--minReadsAtLocus 4 \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-10.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-11.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-12.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-13.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-14.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-15.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-16.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-17.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-18.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-1.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-2.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-3.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-4.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-5.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-6.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-7.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-8.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/ery_30-9.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-10.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-11.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-12.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-13.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-14.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-15.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-16.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-17.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-18.sorted.bam \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-1.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-2.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-3.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-4.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-5.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-6.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-7.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-8.sorted.bam  \
# # 	-I /data3/claudius/Big_Data/BOWTIE2/BAM/par_34-9.sorted.bam  
# #
# #
# #
# #
# # # I am going to create a VCF file with freebayes, which does realignment around indels inherently.
# # # I am going to use the indels it reports in its VCF as argument to the "-known" flag of GATK's
# # # RealignerTargetGenerator in the hope that it will run faster then.
# # cd /data3/claudius/Big_Data/BOWTIE2/BAM
# #
# #
# #
# # # The following command takes the contig ID's and their lengths from a BAM file header. It then
# # # reads ID's and lengths into a Perl hash table. It then opens the *angsd_regions" positive list
# # # of contigs without SE mismappings and SE coverage < Q95. It then prints out those contig ID's
# # # together with their respective lengths. This list is then numerically sorted by contig number
# # # and then turned into BED intervals. BED intervals use 0-based coordinates and the END coordinate
# # # is not included in the interval. Since the contigs are padded with a single N on either end, I
# # # therefore just have to add a "1" to the interval start.
# # samtools view -h ery_30-11.sorted.bam | grep "^@SQ" | sed 's/^.*SN://' | sed 's/LN://' | \
# # 	perl -ane'
# # 	$H{$F[0]}=$F[1];
# # 	END{open($IN, "cat Big_Data_Contigs.noSEgt2.depthLT26.angsd_regions | sed 's/://' |"); 
# # 	while(<$IN>){chomp; die unless exists $H{$_}; print $_, "\t", $H{$_}, "\n";}}
# # 		' | \
# # 	sort -V -k1,1 | gawk '{print $1"\t"1"\t"$2}' > Big_Data_Contigs.noSEgt2.depthLT26.bed
# # 
# # # ==> A test run of an individual BAM file with freebayes resulted in a segmentation fault. This is likely due
# # # to a limitation of freebayes on the number of reference contigs in the BAM file, see issue:
# # # https://github.com/ekg/freebayes/issues/321
# # # freebayes works fine with 200k contigs but not with 250k and the Big Data reference has >500k contigs
# # # ==> I am therefore going to concatenate contigs into pseudo scaffolds.
# 
# # ------ PSEUDO SCAFFOLDS --------
# 
# # cd /data3/claudius/Big_Data/FREEBAYES
# # mkdir contigsToScaffold_ref
# # # ---------------------------------------------
# # # Generate scaffold as workaround for freebayes
# # # ---------------------------------------------
# # # The following command takes the reference contigs (583,312) and concatenates them into long contigs.
# # # It first writes the fasta header of the first scaffold (1000 contigs shall be concatenated into one scaffold).
# # # It then reads in only the sequence lines from the fasta input. It then extracts what should be the remainder of
# # # the SbfI restriction site. It then counts the number of mismatches of this sequence to the expected TGCAGG. If the
# # # the number of mismatches is greater than 2, then skip this contig and move on to the next. Note, those mismatches
# # # should be caused by misassembly, all input reads were filtered for having the correct remainder of the retriction
# # # site by `process_radtags`. The contig sequence is printed out followed by a stretch of N's that makes the contig
# # # 850 bp long. Note, the longest contig is 787 bp long. So this N-padding should prevent read mappings across contigs.
# # # Finally it checks whether the input line is a multiple of 1000 and, if true, ends this scaffold and starts a new one.
# # cat Big_Data_ref.fa | perl -ne'
# # 	BEGIN{print ">Contig_1\n";}
# # 	if($.%2==0){
# # 		($rs) = $_ =~ m/^.(......).*/; 
# # 		$mm = ($rs ^ "TGCAGC") =~ tr/\001-\255//; 
# # 		next if $mm > 2; 
# # 		chomp, $a = length; print; 
# # 		print "N" x (850 - $a); 
# # 		if($.%1000==0){print "\n>Contig_", $./1000+1, "\n"}
# # 		}' > contigsToScaffold_ref/Big_Data_ref_scaffolds.fa
# #
# # Note, the fasta of each scaffold contains the "Contig". This was a choice of naming was misleading. A better fasta header would contain the word "scaffold" here.
# # 
# 
# # # ---------------
# # # bowtie2 index
# # # ---------------
# # cd contigsToScaffold_ref
# # bowtie2-build  Big_Data_ref_scaffolds.fa  Big_Data_ref_scaffolds -f -o 4 -t 12 --threads 23
# # # since huluvu has a lot of memory, I am slightly increasing the offrate and lookup table of the index, 
# # # both should make searching the index faster
# 
# # # ----------------
# # # bowtie2 mapping
# # # ----------------
# # mkdir BAM
# # for SE in ../../data/*fq_1.gz; do
# # 	PE=`echo $SE | sed 's/fq_1/fq_2/'`
# # 	ID=`basename $SE .fq_1.gz` 
# # 	bowtie2 -x Big_Data_ref_scaffolds \
# # 		-1 <(zcat $SE) -2 <(zcat $PE) \
# # 		-q \
# # 		--phred64 \
# # 		--very-sensitive \
# # 		--dpad 10 \
# # 		--gbar 4 \
# # 		--norc \
# # 		--end-to-end \
# # 		--np 10 \
# # 		-I 60 \
# # 		-X 800 \
# # 		--no-contain \
# # 		-t \
# # 		--no-unal \
# # 		--rg-id $ID \
# # 		--rg "SM:$ID" \
# # 		--omit-sec-seq \
# # 		-p 22 \
# # 		| samtools view -uq 1 - \
# # 		| samtools sort -o BAM/$ID.sorted.bam;
# # done
# # # # running bowtie2 in end-to-end mode (i. e. no clipping of query sequences), "very sensitive" specifies a "seed" (kmer) length of 20
# # # # that needs to match exactly, seeds a sampled from the query according to the function f(x) = 1 + 0.5 * sqrt(x), x is the read length
# # # # with "--dpad 10" I am allowing gaps to be up to 10 bp long; "--gbar 4" disallows gaps within 4 bp of either end of the reada; with
# # # # "--norc" bowtie2 will not attempt to align unpaired reads against the revcomp reference strand - I generally don't expect any good alignments
# # # # against the revcomp of the RADome, so this flag somewhat reduces unnecessary searches; "--np 10" sets the penalty for ambiguous 
# # # # characters in an alignment - in the RADome SE RAD tags and PE contigs are separated by up to 10 N's, I don't want alignments across this gap;
# # # # I specify a minimum fragment length of 60 (i. e. allowing some overlap between SE and PE reads) and a maximum fragment length of 800 (longer
# # # # than the longest sequence in the RADome; I am supressing the output of SAM records for reads that failed to align anywhere; if secondary sequences
# # # # are printed (which they shouldn't with this configuration), then omit the seq and qual string; "-p 22" runs bowtie2 on 22 cores; samtools view
# # # # then filters out all SAM records with 0 mapping quality and writes out a compressed BAM file
# # 
# # # print a summary of the distribution of edit distances reported in the NM tag in the SAM output:
# # samtools view par_34-9.sorted.bam  | egrep -o "NM:i:[[:digit:]]+" | sed 's/NM:i://' | Rscript -e 'x=scan("stdin");cat(summary(x),fill=T)'
# # 
# 
# # # ----------------------------------------------------
# # # Mismapping filter on BAM files from pseudo scaffold
# # # ----------------------------------------------------
# # parallel "samtools view -h {} | ./filter_mismapped_reads.pl | samtools view -b - > {.}.MismapFiltered.bam" ::: *sorted.bam
# # # see the filter script 'filter_mismapped_reads.pl' for details of the filter. Briefly, I am filtering SAM records, not whole
# # # contigs, of mismapping SE and PE reads. A mismapping SE read has a mapping position that if reduced by 2 is not a multiple of 850,
# # # the contig length that I set during creation of pseudo scaffolds. A PE read is mismapped if it maps to the SE RAD tag portion of
# # # a contig (as created by 'rbasm', not a scaffold contig). However, I do allow those mappings if it is a proper pair mapping, i. e.
# # # SE and PE reads overlap and the fragment length is >=60, as I set as minimum ISIZE for bowtie2 mapping.
# 
# # # ----------------------------------------------------
# # # DUST filtering of pseudo scaffold reference
# # # ----------------------------------------------------
# # 
# # # find low complexity intervals, takes only 5 minutes
# # dustmasker -in Big_Data_ref_scaffolds.fa -out Big_Data_ref_scaffolds.dust_intervals -window 64 -level 20 -parse_seqids -outfmt acclist
# # 
# # # turn dust intervals to BED format: 
# # sed 's/^>lcl|//' *intervals | gawk 'BEGIN{OFS="\t"}{$3+=1; print}' > Big_Data_ref_scaffolds.dust_intervals.bed
# # 
# # # create index for BEDtools:
# # perl -ne'if($.%2){s/^>//; chomp; print}else{chomp; print "\t", length, "\n"}' Big_Data_ref_scaffolds.fa > Big_Data_ref_scaffolds.idx
# # 
# # # get complement of low complexity intervals, i. e. high complexity intervals:
# # bedtools complement -i Big_Data_ref_scaffolds.dust_intervals.bed -g Big_Data_ref_scaffolds.idx > Big_Data_ref_scaffolds.dust_intervals.compl.bed
# # 
# # # # This still needs to be transformed to ANGSD format:
# # gawk 'BEGIN{OFS="\t"}{$2+=1; print}' Big_Data_ref.dust_intervals.sorted.compl.bed > Big_Data_ref.dust_intervals.sorted.compl.angsd_sites
# # # # ANGSD uses 1-indexing of positions and includes regStart and regEnd in the interval
# 
# 
# # # -------------------------------------------------------
# # # excess coverage filtering of pseudo scaffold reference
# # # -------------------------------------------------------
# # # get 99th percentiles of coverage distributions of SE RAD tags from each ind
# # # note that RAD tags with 0 coverage are not included in the distributions
# # # note that coverage on paired-end contigs is not considered
# # echo "# filename       Q99_cov" > Q99.cov; \
# # for f in *sorted.MismapFiltered.bam; 
# # do 
# # 	echo -n "$f " >> Q99.cov; \
# # 	samtools view -f 64 $f | cut -f 3-4 | uniq -c | gawk '{print $1}' | \
# # 	Rscript -e 'x=scan("stdin",quiet=T);cat(quantile(x,probs=c(0.99)),fill=T);' \
# # 	>> Q99.cov;
# # done
# # 
# # # prints out 'Q99.excess_cov.intervals.bed', the intervals only contain SE RAD tag regions with excess coverage,
# # # a region has excess coverage if in any of the individual BAM files it has coverage above the 99th percentile
# # # of the individuals coverage distribution
# # ./print_Q99_excessCov_intervals.pl
# # # this takes Q99_cov as input and spits out Q99.excess_cov.intervals.bed
# # 
# # # Now I need to complement those regions, as done above for the dust intervals
# # sort -V -k1,1 -k2,2 Q99.excess_cov.intervals.bed > Q99.excess_cov.intervals.sorted.bed 
# # bedtools complement -i Q99.excess_cov.intervals.sorted.bed -g ../Big_Data_ref_scaffolds.idx > Q99.excess_cov.intervals.sorted.compl.bed
# # 
# #
# #
# # # -----------
# # # FREEBAYES
# # # -----------
# # # run freebayes in parallel
# # # I am only going to limit freebayes to the no-excess-coverage BED intervals, not the dust intervals, since I want freebayes to call indels
# # # at low complexity sites. I am parallelising freebayes with GNU parallel as shown in the freebayes-parallel script in the freebayes repo.
# # # I am going to execute freebayes with different sets of BED intervals as target regions. So, I first need to split the no-excess-coverage BED
# # # interval file:
# # mkdir split_BEDs
# # split -l 60 Q99.excess_cov.intervals.sorted.compl.bed split_BEDs/
# # # this will output arbitrarily named files, each containing 60 lines from the BED file
# # 
# # # freebayes needs BAM indeces to work:
# # ls *MismapFiltered.bam | parallel samtools index {}
# # 
# # # run different processes of freebayes with different target regions in parallel:
# # ls split_BEDs/* | \
# # 	/usr/bin/time -vo freebayes_ParEry.time \
# # 	parallel -k -j 12 "freebayes -L bamfile.list -f ../Big_Data_ref_scaffolds.fa -t {} --populations pop.list -T 0.01 -n 6 -m 2 -F 0.2 --min-coverage 5 -V -j" \
# # 	| vcffirstheader | vcfstreamsort | vcfuniq | bgzip > freebayes_ParEry.vcf.gz
# # # I am running freebayes with all 36 BAM files in the bamfile.list. I also inform freebayes that the individuals are from two groups and that therefore HWE should
# # # only be assumed for each group separately, not for both groups together. I am also setting a prior expected nucleotide heterozygosity (pi) of 0.01, limiting the
# # # the number of alleles to evaluate at site to 6 (in the hope that it does not refer to SNP's, which of course have only 4 alleles), require a min mapping quality 
# # # of 2, require a min alternate allele fraction of 0.2 (default) to consider a site and a min coverage of 5 (I am assuming this means coverage across ind, not per ind).
# # # I also have to turn off some prior expectations (-V) about the read data that are violated in RAD data. Finally, I am specifying to include mapping quality into
# # # the calculation of likelihoods.
# # # This run took only 21min.
# # 
# # # --------------------------
# # # create indelOnly VCF file  
# # # --------------------------
# # # VCF file needs to BGZIPed and indexed with tabix
# # zcat freebayes_ParEry.vcf.gz | vcfindels | vcffilter -sf "DP > 10" | vcffilter -sf "QUAL > 20" | bgzip > freebayes_ParEry_indelsOnly_QualGt20.vcf.gz
# # # the filter scripts are from vcflib, I filter out sites with 10 or less coverage and 20 or less Phred scaled prob. of being polymorphic.
# # 
# # # ----------------------------
# # # indel realignment with glia
# # # ----------------------------
# # /usr/bin/time -vo glia_realign.time \
# # 	bamtools merge -list bamfile.list | \
# # 	glia -f ../Big_Data_ref_scaffolds.fa -R -v freebayes_ParEry_indelsOnly_QualGt20.vcf.gz -w 800 -Q 200 -C 2 | 
# # 	bamtools split -tag RG -tagPrefix "" -stub "realigned" 
# # # the indel realignment took 20 minutes
# # # with upper command line I am first merging all individual BAM files, then realigning them mostly around the indels in the VCF file
# # # and then splitting the BAM stream again by read group. I am providing glia with a merged BAM stream because I hope that glia doesn't just
# # # do indel realignment on a read-by-read basis but rather does something similar to a local assembly of reads. For the glia command I specify a window size
# # # of 800 bp, forcing realignment of reads if the sum of their qualities of mismatched bases is greater than 200 or if the number of gaps is greater than 2.
# # # Glia realignment has minimal to no effect on the alignment of reads. I used commands like:
# # # diff <(samtools view ../par_34-9.sorted.MismapFiltered.bam Contig_30:) <(samtools view realigned_par_34-9.bam Contig_30:) | less -S
# # # to look at the difference between BAM alignments before and after glia realignment. I could not see any change in mapping position and no change in the number
# # # or size of deletion as reported in the CIGAR string. The only difference I could see was the introduction of the CIGAR symbol 'X' for a mismatch. Usually 
# # # mapping programmes output an 'M' which stands for an *alignment* match, but could be a match or a mismatch.
# # # ==> glia realignment had no significant effect
# # 
# # # ---------------
# # # BED intervals
# # # ---------------
# # # I want to create an intervals BED file that contains all intervals that passed the Q99 coverage filter from above and NOT the low complexity intervals
# # # detected with dustmaker (see above):
# # bedtools subtract -a ../Q99.excess_cov.intervals.sorted.compl.bed -b ../../Big_Data_ref_scaffolds.dust_intervals.bed > Big_Data_ref_scaffolds_noExCov_noDUST.intervals.bed
# # 
# # # the total length of the new interals:
# # gawk 'BEGIN{l=0}{l+=($3-$2)}END{print l}' Big_Data_ref_scaffolds_noExCov_noDUST.intervals.bed
# # # this prints: 89114979 , i. e. ~89Mb
# # 
# # # -----------------------------------------------------
# # # run freebayes in parallel on all realigned BAM files:
# # # -----------------------------------------------------
# # ls split_BEDs/* | \
# # 	/usr/bin/time -vo freebayes_ParEry.time \
# # 	parallel -k -j 12 "freebayes -L bamfile.list -f ../Big_Data_ref_scaffolds.fa -t {} --populations pop.list -T 0.01 -n 6 -m 2 -F 0.2 --min-coverage 5 -V -j" \
# # 	| vcffirstheader | vcfstreamsort | vcfuniq | bgzip > freebayes_ParEry.vcf.gz
# # # remember to use 'bgzip' instead 'gzip' to compress the VCF file !!!
# # # this command takes 50 min to complete
# # # Freebayes by default does not report non-variable sites, but if I want to use the VCF for filtering and subsequent conversion into a "sites.keep" file for ANGSD,
# # # I need to have all sites in the BED intervals in the VCF output.
# # 
# # ls split_BEDs/* | \
# # 	/usr/bin/time -vo freebayes_ParEry_inclMono_realigned.time \
# # 	parallel -k -j 20 "freebayes -L bamfile.list -f ../../Big_Data_ref_scaffolds.fa -t {} \
# # 	--populations pop.list -T 0.01 -n 6 -m 5 -F 0.2 --min-coverage 10 -V -j --report-monomorphic" | \
# # 	vcffirstheader | vcfstreamsort | vcfuniq | bgzip \
# # 	> freebayes_ParEry_inclMono_realigned.vcf.gz & 
# # # # Note, running frebayes with --report-monomorphic takes much longer, almost 4h ! Most of the time freebayes spends on the first 10,000 intervals. I had split
# # # # the "noExCov_noDUST" intervals into files of 10,000 intervals each. So only one freebayes process was working on the first split file. I have therefore split
# # # # that file again on every 100th interval to allow parallel processing of these first intervals which seemingly got much higher coverage than others.
# # # # When freebayes or samtools mpileup are parallelised like this, it is important that the BED intervals are fed into GNU parallel in the right order, i. e. chromosome and position
# # # # sorted. That means that after splitting a large sorted BED file with Unix "split", the split BED file with the first few intervals needs to have a file name that comes first in lexical
# # # # order, as output by Unix "ls" by default. If this is not taken care of, then it is not guaranteed that the final VCF file will be position sorted, which is necessary
# # # # for vcflib's "vcffilter" to work. So when splitting a second time, e. g. the file "aa" produced by the previous split coomand, then use "aa" as the prefix for the output file
# # # # names of the second split command, which will produce files with names: aaaa, aaab, aaac ...
# # 
# # # --------
# # # ANGSD
# # # --------
# # 
# # # turn the BED file into ANGSD region format
# # cat Big_Data_ref_scaffolds_noExCov_noDUST.intervals.bed | \
# # 	gawk '{print $1"\t"$2+1"\t"$3}' | \
# # 	gawk '{print $1":"$2"-"$3}' \
# # 	>Big_Data_ref_scaffolds_noExCov_noDUST.intervals.angsd_regions
# # # I tried to run angsd with the more than 1.2M regions, but it takes too long to read in data for those regions. I am instead going to let ANGSD
# # # read in all BAM records but limit its analyses to certain sites in the hope that this will be faster.
# # 
# # # Since freebayes "--report-monomorphic" produces much fewer VCF records as there are sites in the BED intervals,
# # # I am now trying to create a genomic VCF with freebayes. I am also specifying a haplotype length of 1, in order
# # # to suppress the output of haplotypes.
# # ls split_BEDs/* | \
# # 	/usr/bin/time -vo freebayes_ParEry_inclMono_realigned.time \
# # 	parallel -k -j 12 "freebayes -L bamfile.list -f ../../Big_Data_ref_scaffolds.fa -t {} \
# # 	--populations pop.list --gvcf --gvcf-chunk 1 -T 0.01 -n 6 --haplotype-length 1 -m 5 -F 0.2 -C 2 -V -j" | \
# # 	vcffirstheader | vcfstreamsort | vcfuniq | bgzip \
# # 	> freebayes_ParEry_realigned.gvcf.gz & 
# # # This command takes 1h to complete. However, non-variable reference sites are output in small blocks, almost all larger than 1 base long,
# # # as opposed what I thought "--gvcf-chunk 1" would do.
# 
# # # -----------------
# # # SAMTOOLS MPILEUP
# # # -----------------
# # # Samtools mpileup seems to output a VCF record for every position in the reference sequence, i. e. a gVCF file.
# # cd /data3/claudius/Big_Data/FREEBAYES/contigsToScaffold_ref/BAM/glia_realigned
# # ls split_BEDs/* | \
# # 	/usr/bin/time -vo samtools.time \
# # 	parallel -k -j 12 \
# # 	"samtools mpileup -l {} -s -Q 20 -A -t DP -d 5000 -uv -f ../../Big_Data_ref_scaffolds.fa realigned_*bam" \
# # 	| vcffirstheader | vcfstreamsort | vcfuniq \
# # 	| bcftools call -m -P 0.01 -O z -o samtools_mpileup_ParEry_realigned.vcf.gz
# # # this takes 1h 38min to complete
# # # this gives all glia realigned BAM files to samtools mpileup
# # # "-s" outputs root mean square mapping quality at the MQ tag in the INFO column
# # # to bcftools I am giving a prior mutation rate of 0.01 and with -m I am turning on joint variant calling across all individuals (I hope, not clearly documented)
# # # there seems to be a strong bias towards deletions from the reference
# # # About 99.4% of sites in the BED file are included in the VCF that the upper command line produces.
# # # NOTE: samtools mpileup also calls outside the intervals that are fed into it in the above commands !!! So DUST and excess coverage filtering has not effect.
# # #       So samtools/bcftools is caling indels in repetitive sequences.
# # 
# # # -------- 
# # # GATK
# # # --------
# # # Since samtools/bcftools seems to be buggy (see above), I am going to try out GATK for creating gVCF that I could use for filtering.
# # #
# # # I am first realigning reads around known indels infered with freebayes (see above). 
# # # 1) creating target intervals:
# # cd /data3/claudius/Big_Data/FREEBAYES/contigsToScaffold_ref/BAM
# # java -jar $gatk \
# # 	-T RealignerTargetCreator \
# # 	-R Big_Data_ref_scaffolds.fa \
# # 	-known freebayes_ParEry_indelsOnly_QualGt20.vcf \
# # 	-I ParEry.sorted.MismapFiltered.merged.bam \
# # 	-o ParEry.sorted.MismapFiltered.merged.RealignerTargetIntervals
# # # Note, that stupidly GATK needs the ref fasta, .fai and .dict files in the same directory as where the command is executed from. The files can also not be compressed.
# # # I am providing a VCF file with indels to force GATK to create intervals around those. I am also providing the merged read data from all inds (-I). Intervals are also
# # # infered from CIGAR strings in the BAM file:
# # # "The tool adds indel sites present in the known indels file and indel sites in the alignment CIGAR strings to the targets. 
# # # Additionally, the tool considers the presence of mismatches and soft-clips, and adds regions that pass a concentration threshold to the target intervals."
# # # (from: http://gatkforums.broadinstitute.org/gatk/discussion/7156/howto-perform-local-realignment-around-indels#latest )
# # # This creates 123,437 intervals in ~25 min
# # #
# # # 2) indel realignment:
# # /usr/bin/time -vo IndelRealinger.time \
# # 	java -Xmx90G -jar $gatk \
# # 	-T IndelRealigner \
# # 	-R Big_Data_ref_scaffolds.fa \
# # 	-targetIntervals ParEry.sorted.MismapFiltered.merged.RealignerTarget.intervals \
# # 	-I ParEry.sorted.MismapFiltered.merged.bam \
# # 	-model USE_READS \
# # 	-maxReads 3000 \
# # 	-maxConsensuses 5 \
# # 	-o ParEry.sorted.MismapFiltered.merged.GATK_IndelRealigned.bam
# # # Note that I had to rename the file containing the target intervals, because stupidly GATK IndelRealigner requires either the file ending ".intervals" or ".list". 
# # # The USE_READS consensus model constructs alternative alignments from the reference sequence from indels in reads spanning the site.
# # # I am ignoring sites with coverage above 3000 and let GATK only consider up to 5 possible alternative consensuses to align reads to. Both should speed the process,
# # # IndelRealigner cannot be parallelised.
# # # The resulting BAM file contains unchanged and modified BAM reacords. Realigned reads have BAM records with a new BAM TAG: OC, that contains the old CIGAR string.
# # # The indel realignment is finished after ~32 min.
# # # I checked for differences between BAM records before and after GATK realignment with commands like:
# # diff <(samtools view ../par_34-9.sorted.MismapFiltered.bam Contig_99:) <(samtools view par_34-9_GATKrealigned.bam Contig_99:) | less -S
# # # GATK's IndelRealigner should add a new SAM tag for a read it realigned containing the old CIGAR string: 
# # samtools view par_34-1_GATKrealigned.bam | grep "OC:" | less -S
# # # This finds nothing. There are no realignments what so ever !!!
# # #
# # # 3) split realigned BAM file by read group (i. e. ind):
# # samtools split -f 'GATK_realigned/%!_GATKrealigned.%.' ParEry.sorted.MismapFiltered.merged.bam
# # 
# # # ---------------
# # # HaplotypeCaller
# # # ---------------
# # #
# # cd /data3/claudius/Big_Data/FREEBAYES/contigsToScaffold_ref/BAM
# # # 1) convert BED format of reference intervals to 1-based closed intervals needed for GATK:
# # gawk '{print $1":"$2+1"-"$3}' Big_Data_ref_scaffolds_noExCov_noDUST.intervals.bed > Big_Data_ref_scaffolds_noExCov_noDUST.intervals
# # #
# # # 2) run GATK HaplotypeCaller on all individual BAM files in paralel:
# # ls GATK_realigned/*bam | \
# # 	/usr/bin/time -vo GATK_HaplotypeCaller.time \
# # 	parallel -j 6 "java -Xmx10G -jar $gatk -T HaplotypeCaller -R Big_Data_ref_scaffolds.fa -I {} \
# # 	-L Big_Data_ref_scaffolds_noExCov_noDUST.intervals -ERC BP_RESOLUTION -mmq 5 -mbq 10 -minReadsPerAlignStart 0 -inputPrior 0.33 -inputPrior 0.33  -o {.}.g.vcf.gz"
# # # This command took 1h 40 min to finish.
# # # Thise creates a gVCF for every individual separately. 
# # # Note: for GATK the naming of the output file is very important. The output file needs to end either in .vcf or .vcf.gz, not .gvcf. If the output file name ends with
# # # .vcf.gz, then GATK will automatically bgzip it and create a tabix index. Bgzipping myself (with process substitution) and creating an index separately with tabix leads
# # # to errors with downstream programmes. This might be due to different kinds/version of index created or compression used.
# # # I am running HAplotypeCaller in EmitRefConfidence (ERC) mode, that emits a genomic VCF, i. e. a VCF with a record for every position in the reference. 
# # # I am requiring a min mapQ of 5 and I specify a flat prior on the alternate allele counts 0, 1 and 2. Note, this makes much more sense than specifying
# # # a heterozygosity prior, i. e. theta, since with a flat prior I am not biasing variant/genotype calling for the reference sequence. That is, the ref
# # # sequence the way I created it with rainbow and Vsearch does not necessarily always contain the ancestral allele at every position.
# # #
# # 3) Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file:
# # # Note: this takes too long. Skip to 4)
# # cd /data3/claudius/Big_Data/FREEBAYES/contigsToScaffold_ref/BAM/GATK_realigned
# # # variant(){ ls *vcf.gz | perl -ne'chomp; print "--variant $_ "' }
# # java -jar -Xmx50G $gatk -T CombineGVCFs -R ../Big_Data_ref_scaffolds.fa \
# # 	--variant ery_30-10_GATKrealigned.g.vcf.gz --variant ery_30-11_GATKrealigned.g.vcf.gz --variant ery_30-12_GATKrealigned.g.vcf.gz --variant ery_30-13_GATKrealigned.g.vcf.gz --variant ery_30-14_GATKrealigned.g.vcf.gz --variant ery_30-15_GATKrealigned.g.vcf.gz --variant ery_30-16_GATKrealigned.g.vcf.gz --variant ery_30-17_GATKrealigned.g.vcf.gz --variant ery_30-18_GATKrealigned.g.vcf.gz --variant ery_30-1_GATKrealigned.g.vcf.gz --variant ery_30-2_GATKrealigned.g.vcf.gz --variant ery_30-3_GATKrealigned.g.vcf.gz --variant ery_30-4_GATKrealigned.g.vcf.gz --variant ery_30-5_GATKrealigned.g.vcf.gz --variant ery_30-6_GATKrealigned.g.vcf.gz --variant ery_30-7_GATKrealigned.g.vcf.gz --variant ery_30-8_GATKrealigned.g.vcf.gz --variant ery_30-9_GATKrealigned.g.vcf.gz --variant par_34-10_GATKrealigned.g.vcf.gz --variant par_34-11_GATKrealigned.g.vcf.gz --variant par_34-12_GATKrealigned.g.vcf.gz --variant par_34-13_GATKrealigned.g.vcf.gz --variant par_34-14_GATKrealigned.g.vcf.gz --variant par_34-15_GATKrealigned.g.vcf.gz --variant par_34-16_GATKrealigned.g.vcf.gz --variant par_34-17_GATKrealigned.g.vcf.gz --variant par_34-18_GATKrealigned.g.vcf.gz --variant par_34-1_GATKrealigned.g.vcf.gz --variant par_34-2_GATKrealigned.g.vcf.gz --variant par_34-3_GATKrealigned.g.vcf.gz --variant par_34-4_GATKrealigned.g.vcf.gz --variant par_34-5_GATKrealigned.g.vcf.gz --variant par_34-6_GATKrealigned.g.vcf.gz --variant par_34-7_GATKrealigned.g.vcf.gz --variant par_34-8_GATKrealigned.g.vcf.gz --variant par_34-9_GATKrealigned.g.vcf.gz \
# # 	-o ParEry_combined.g.vcf.gz
# # 
# # # 4) Perform joint genotyping on gVCF files produced by HaplotypeCaller:
# # # The following let to error messages:
# #	# cd /data3/claudius/Big_Data/FREEBAYES/contigsToScaffold_ref/BAM/GATK_realigned
# #	# variant(){ 
# #	# 	ls *vcf.gz | perl -ne'chomp; print "-V $_ "' 
# #	# }
# #	# export -f variant
# #	# prior(){ 
# #	# 	perl -e'$p=1/(72+1); print "-inputPrior $p " x 72' 
# #	# }
# #	# export -f prior
# #	# /usr/bin/time -vo GenotypeGVCFs.time java -jar $gatk \
# #	# 	-T GenotypeGVCFs \
# #	# 	-R ../Big_Data_ref_scaffolds.fa \
# #	# 	$(variant) \
# #	# 	-nt 6 \
# #	# 	-stand_emit_conf 0 \
# #	# 	-stand_call_conf 20 \
# #	# 	-allSites \
# #	# 	$(prior) \
# #	# 	-o ParEry.vcf.gz \
# #	# 	2> GenotypeGVCFs.log
# #	# # This produces the following error message:
# #	# # ERROR MESSAGE: Invalid command line: Argument inputPrior has a bad value: Invalid length of inputPrior vector: vector length must be equal to # samples +1
# # # The following command works and runs in parallel mode:
# # cd /data3/claudius/Big_Data/FREEBAYES/contigsToScaffold_ref/BAM/GATK_realigned
# # /usr/bin/time -vo GenotypeGVCFs.time java -jar $gatk \
# # 	-T GenotypeGVCFs \
# # 	-R ../Big_Data_ref_scaffolds.fa \
# # 	-nt 22 \
# #  	-stand_call_conf 20 \
# #  	-stand_emit_conf 0 \
# #  	-allSites \
# # 	-V ery_30-10_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-11_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-12_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-13_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-14_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-15_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-16_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-17_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-18_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-1_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-2_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-3_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-4_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-5_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-6_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-7_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-8_GATKrealigned.g.vcf.gz \
# # 	-V ery_30-9_GATKrealigned.g.vcf.gz \
# # 	-V par_34-10_GATKrealigned.g.vcf.gz \
# # 	-V par_34-11_GATKrealigned.g.vcf.gz \
# # 	-V par_34-12_GATKrealigned.g.vcf.gz \
# # 	-V par_34-13_GATKrealigned.g.vcf.gz \
# # 	-V par_34-14_GATKrealigned.g.vcf.gz \
# # 	-V par_34-15_GATKrealigned.g.vcf.gz \
# # 	-V par_34-16_GATKrealigned.g.vcf.gz \
# # 	-V par_34-17_GATKrealigned.g.vcf.gz \
# # 	-V par_34-18_GATKrealigned.g.vcf.gz \
# # 	-V par_34-1_GATKrealigned.g.vcf.gz \
# # 	-V par_34-2_GATKrealigned.g.vcf.gz \
# # 	-V par_34-3_GATKrealigned.g.vcf.gz \
# # 	-V par_34-4_GATKrealigned.g.vcf.gz \
# # 	-V par_34-5_GATKrealigned.g.vcf.gz \
# # 	-V par_34-6_GATKrealigned.g.vcf.gz \
# # 	-V par_34-7_GATKrealigned.g.vcf.gz \
# # 	-V par_34-8_GATKrealigned.g.vcf.gz \
# # 	-V par_34-9_GATKrealigned.g.vcf.gz \
# # 	-o ParEry.vcf.gz \
# # 	2> GenotypeGVCFs.log
# # # Questions: The parameter --input_prior can be applied both at the HaplotypeCaller stage and at the GenotypeGVCFs stage. 
# # # I'm wondering if it makes any difference at what stage this parameter is applied?
# # # Answer: It will only make a difference in HaplotypeCaller, not GenotypeGVCFs. The option is available for GenotypeGVCFs 
# # # because of how the genotyping tools inherit arguments from shared code classes. But, again, the argument does nothing in GenotypeGVCFs.
# # # Taken from: http://gatkforums.broadinstitute.org/gatk/discussion/7370/input-prior-applied-at-haplotypecaller-vs-genotypegvcfs-stage
# # # ==> So no 'inputPrior' for the tool GenotypeVCFs, only for Haplotype caller !
# # 
# # # ------------------------------
# # # Analysis of the ParEry.vcf.log
# # # ------------------------------
# # # + ALT allele count (AC) is always 0 on Contig_21 (scaffold), i. e. no alternate allele is called at all !
# # # + Contig_11	1762	.	A	G,<NON_REF>	139.96	.	AC=0,0;AF=0.00,0.00;AN=2;DP=3;MLEAC=0,0;MLEAF=0.00,0.00;MQ=9.82	GT:AD:DP:RGQ	./.:0,0,0:0	./.:0,0,0:./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	0/0:1,0,0:1:3	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0
# # # 	- total depth allegedly 3, but only one called genotype, it reports one read, for ref allele, still the ALT allele G has been called with 139.96 ?!
# # # + Contig_11	1913	.	C	T,<NON_REF>	64.25	.	AC=0,0;AF=0.00,0.00;AN=16;DP=111;ExcessHet=3.01;InbreedingCoeff=0.6674;MLEAC=0,0;MLEAF=0.00,0.00;MQ=6.87	GT:AD:DP:RGQ	./.:7,0,0:7	0/0:8,0,0:8:3	0/0:8,0,0:8:6	./.:8,0,0:8	./.:3,0,0:3	./.:2,0,0:2	./.:2,0,0:2	0/0:9,0,0:9:24	0/0:7,0,0:7:3	./.:2,0,0:2	./.:3,0,0:3	./.:4,0,0:4	./.:10,0,0:10	./.:0,0,0:0	0/0:2,0,0:2:6	0/0:2,0,0:2:3	0/0:2,0,0:2:3	./.:7,0,0:7	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	./.:4,0,0:4	./.:0,0,0:0	./.:1,0,0:1	./.:1,0,0:1	./.:0,0,0:0	./.:0,0,0:0	./.:0,0,0:0	0/0:1,0,0:1:3	./.:0,0,0:0	./.:1,0,0:1	./.:0,0,0:0	./.:0,0,0:0
# # # 	- take the first genotype/sample column, it has 7 reads for the ref allele, but no genotype called, the whole locus has 111 coverage, although that does
# # # 	not match with the sum of allele depth reported in the genotype columns, so I think the total depth must be including reads that have been filtered
# # # + Contig_11	4437	.	C	T,<NON_REF>	360.46	.	AC=5,0;AF=0.100,0.00;AN=50;BaseQRankSum=-1.534e+00;ClippingRankSum=0.00;DP=202;ExcessHet=0.0154;FS=0.000;InbreedingCoeff=0.6387;MLEAC=4,0;MLEAF=0.080,0.00;MQ=25.36;MQRankSum=-1.534e+00;QD=27.73;ReadPosRankSum=-1.570e-01;SOR=1.508	GT:AD:DP:GQ:PL	./.:0,0,0:0:.:0,0,0,0,0,0	0/0:14,0,0:14:42:0,42,609,42,609,609	0/0:10,0,0:10:30:0,30,428,30,428,428	0/0:3,0,0:3:9:0,9,133,9,133,133	./.:0,0,0:0:.:0,0,0,0,0,0	./.:0,0,0:0:.:0,0,0,0,0,0	0/0:13,0,0:13:39:0,39,555,39,555,555	0/0:7,0,0:7:21:0,21,312,21,312,312	0/0:1,0,0:1:3:0,3,36,3,36,36	0/0:10,0,0:10:30:0,30,444,30,444,444	0/0:14,0,0:14:42:0,42,603,42,603,603	0/0:13,0,0:13:33:0,33,495,33,495,495	0/0:2,0,0:2:6:0,6,89,6,89,89	./.:0,0,0:0:.:0,0,0,0,0,0	0/0:2,0,0:2:6:0,6,90,6,90,90	0/0:10,0,0:10:30:0,30,425,30,425,425	./.:0,0,0:0:.:0,0,0,0,0,0	./.:0,0,0:0:.:0,0,0,0,0,0	0/0:7,0,0:7:21:0,21,277,21,277,277	0/0:6,0,0:6:18:0,18,236,18,236,236	1/1:0,6,0:6:18:168,18,0,168,18,168	0/0:11,0,0:11:12:0,12,180,12,180,180	./.:0,0,0:0:.:0,0,0,0,0,0	0/1:7,1,0:8:6:6,0,243,27,246,273	0/0:8,0,0:8:24:0,24,350,24,350,350	0/0:8,0,0:8:24:0,24,344,24,344,344	0/0:2,0,0:2:6:0,6,88,6,88,88	0/0:9,0,0:9:27:0,27,376,27,376,376	./.:0,0,0:0:.:0,0,0,0,0,0	0/0:10,0,0:10:30:0,30,443,30,443,443	./.:7,0,0:7:.:0,0,0,0,0,0	1/1:0,7,0:7:21:248,21,0,248,21,248	./.:0,0,0:0:.:0,0,0,0,0,0	./.:0,0,0:0:.:0,0,0,0,0,0	0/0:9,0,0:9:27:0,27,378,27,378,378	0/0:2,0,0:2:6:0,6,71,6,71,71
# # # + there are 89,115,735 records in the VCF, which is only a tiny bit more than in the BED intervals.
# # # ==> I cannot use this VCF for filtering
# # 
# # # ---------
# # # Coverage
# # # ---------
# # samtools depth -aa -b Big_Data_ref_scaffolds_noExCov_noDUST.intervals.bed -Q 5 *MismapFiltered.bam | bgzip > ParEry.AllSites.depth.gz
# # # This prints out a tab delimited file with Contig<TAB>POS<TAB>Cov<TAB>Cov<TAB>Cov<TAB>... so here 36 coverage columns.
# # # only with '-aa' does samtools output all positions from the BED intervals
# # # I am only counting reads with mapping quality greater than 5.
# 
# 
# ## ---------------- GOING BACK TO INDIVIDUAL CONTIGS -------------------
# 
# # I have decided not to use a VCF file for filtering. Since I do not need to create a VCF file with FreeBayes, I can do filtering on
# # the unscaffolded reference.
# cd /data3/claudius/Big_Data/BOWTIE2/BAM
# 
# # # --------------------
# # # Recap of the filters
# # # --------------------
# # #
# # # MISMAPPING
# # #
# # # line 250 shows the command line that created a negative list of contigs that had SE reads mapping beyond the position 2:
# # # exclude_contigs.SEposgt2
# # # I have renamed that file to Big_Data_contigs.SEposgt2
# # #
# # #
# # # LOW COMPLEXITY
# # #
# # # line 271 onwards shows how I detected low complexity sequences in the Big Data reference with dustmasker and then turned those LC intervals
# # # into BED format
# # #
# # #
# # # EXCESS COVERAGE
# # #
# # # line 309 - 312 shows the command lines I used to determine the coverage of SE reads for each contig for each individual separately. They
# # # created a .SE_depths file for each ind BAM file. However, in the following I am recalculating SE coverage counts for each individual and
# # # streaming those into R to determine the 99th percentile of the coverage distribution for each individual:
# # cd /data3/claudius/Big_Data/BOWTIE2/BAM
# # echo "# filename       Q99_cov" > Q99.cov; \
# # for f in *sorted.bam; 
# # do 
# # 	echo -n "$f " >> Q99.cov; \
# # 	samtools view -f 64 $f | cut -f 3 | uniq -c | gawk '{print $1}' | \
# # 	Rscript -e 'x=scan("stdin",quiet=T);cat(quantile(x,probs=c(0.99)),fill=T);' \
# # 	>> Q99.cov;
# # done
# # #
# # # Next I want to use these 99 percentiles of the individual coverage distributions in order to detect contigs that have SE read coverage
# # # above the individual Q99 in any individual. The following script does that.
# # print_Q99_exCov_contigs.pl
# # # This takes Q99.cov as input and spits out a list of contigs to exclude due to excessive coverage:
# # # Big_Data_Contigs.gtQ99Cov
# # # A contig has excess coverage if in any of the individual BAM files it has coverage above the 99th percentile
# # # of the individuals coverage distribution. 
# #
# # # I now need to create a BED file with intervals that pass all three filters above, i. e. a positive list.
# # # With the following command I am creating a BED file with intervals comprising the whole Big Data reference, 
# # # except for the initial and terminal N of each contig:
# # samtools view -h ery_30-10.sorted.bam | \
# # 	grep "^@SQ" | \
# # 	sed 's/.*SN://' | \
# # 	sed 's/LN://' | \
# # 	gawk '{print $1"\t"1"\t"$2-1}' | \
# # 	sort -Vk1,1 \
# # 	> Big_Data_ref.fa.bed
# # # This takes the header of a BAM file, which contains ref contigs names and their lengths. I then let the start of each interval spanning
# # # a contig be 1 (i. e. 2nd pos in 0-based coord) and the end be congig-length-1. So the terminal N is excluded from the BED interval.
# # # There are 583,312 such intervals. This is as expected the number of fasta records in Big_Data_ref.fa.
# # #
# # # With the command on line 255, I previously created a positive list of contigs with read mappings. So, I have three files with contig names
# # # to keep or delete from Big_Data_ref.fa.bed:
# # # Big_Data_Contigs_with_mappings.list
# # # Big_Data_Contigs.SEposgt2
# # # Big_Data_Contigs.gtQ99Cov
# # # 
# # # The Perl script subtract.pl takes these four files and outputs a new BED file containing intervals that passed the mismapping and excess
# # # coverage filters:
# # ./subtract.pl | sort -Vk1,1  > Big_Data_Contigs.noSEgt2.nogtQ99Cov.bed
# # #
# # # 583312 Big_Data_ref.fa.bed
# # # 575961 Big_Data_Contigs.noSEgt2.nogtQ99Cov.bed
# # #   7351 filtered out
# # # So far I have filtered 7351 contigs from the original reference. 
# #
# # # Now I need to subtract the DUST intervals from this BED file:
# # bedtools subtract -a Big_Data_Contigs.noSEgt2.nogtQ99Cov.bed -b ../Big_Data_ref.dust_intervals.sorted.bed > Big_Data_Contigs.noSEgt2.nogtQ99Cov.noDUST.bed
# #
# # The total length of the new filtered intervals is:
# # gawk '{sum+=($3-$2)}END{print sum}' Big_Data_Contigs.noSEgt2.nogtQ99Cov.noDUST.bed
# # 88,097,045
# # 
# #
# # EVEN COVERAGE 
# #
# # I want to filter individual sites for sufficiently even coverage across individuals. 
# samtools depth -aa -b Big_Data_Contigs.noSEgt2.nogtQ99Cov.noDUST.bed -Q 5 *sorted.bam | bgzip > ParEry.noSEgt2.nogtQ99Cov.noDUST.depth.gz
# # This will print out a line for each position within the BED intervals. First column is contig name, second is position, the remaining
# # 36 columns contain individual coverage for the site, from reads with mapQ > 5.
# # I have checked that to number of sites in the .depth.gz file is the same as the length of the intervals in the BED file.
# 
# # Finally I am running a filter over this depth table, only keeping sites with 3 or more coverage in at least 15 of the individuals:
# ./even.depth.pl 3 15 ParEry.noSEgt2.nogtQ99Cov.noDUST.depth.gz
# # This should produce the file called:
# # ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.depth.gz
# # This file just contains two columns: contig name <TAB>  position
# # The positions are 1-based. There are 2,038,050 sites in this file.
# # From the 2nd till 7th position of each contig, positions cannot be polymorphic (TGCAGG) which would distort SFS.
# # So I am going to filter out those positions:
# zcat ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.depth.gz | sort -Vk1,1 -k2,2 | gawk '$2 !~ /^[234567]$/' > ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites
# # 
# # I want to know how many RAD contigs do these sites belong to:
# cut -f 1 ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites | uniq | wc -l
# # 34,967
# # Before HWE filtering, how many sites have passed filtering so far:
# wc -l ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites
# # 1,829,329
# #
# #
# # HWE FILTERING
# #
# cd /data3/claudius/Big_Data/ANGSD
# # index sites to keep
# angsd sites index ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites
# # ANGSD requires the specification of a region file in order to work. The format is as described here:
# # http://popgen.dk/angsd/index.php/Input#BAM.2FCRAM
# # This should be same format as for samtools regions.
# #
# # To speed up file reading with ANGSD, I am creating new BAM files, that only contain SAM records for
# # reads that mapped to one of the filtered contigs:
# cd /data3/claudius/Big_Data/ANGSD/Data
# ./slim.BAMs.pl
# cd ..
# # during the BAM conversion of the output of this programme, SAMtools emits some warnings of reads whose mate maps to a contig
# # that is not represented in the header. I think, this should not be a problem.
# #
# # Parallelisation within ANGSD does not work properly. So I am going to apply the same trick as used to
# # parallelise Freebayes. First I split the regions file into small files, each containing only 500 contig names:
# split -l 500 ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.rf split_rf/
# # Now I can spawn a separate ANGSD process for each of these smaller regions files:
# ls split_rf/* | \
# 	parallel -j 12 "angsd -HWE_pval_F 1 -b slim.bamfile.list \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites \
# 	-rf {} -out Results/{/} -doMaf 1 -doMajorMinor 1 -GL 1 -P 2 -minQ 20 \
# 	-only_proper_pairs 0 -minMapQ 5"
# # This takes less than 2 min to complete. I have checked that the combined number of lines in the *hweF.gz outfiles
# # is the same as in the *sites file. Unfortunately, the lowest F value reported is 0. So I cannot use this for paralog
# # filtering.
# #
# # ANGSD has another command for HWE testing, which claims to also test for negative deviations of F:
# ls split_rf/* | parallel -j 12 "angsd -doSnpStat 1 -HWE_pval 1 -b slim.bamfile.list \
# 	-only_proper_pairs 0 -minMapQ 5 -sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites 
# 	-rf {} -out Results/{/} -doMaf 1 -doMajorMinor 1 -GL 1 -P 2"
# # However, the output of this command cannot be right. F is always 0.666667, LRT is negative, pval is always 1 and the estimated
# # allele frequencies are much too low. So this function also cannot be used for HWE filtering.
# #
# # I have found another programme, vt, from Goncalo Abecasis' group, that should also calculate HWE and inbreeding statistics
# # from genotype likelihoods. This programme takes genotype likelihoods in a VCF file, i. e. does not calculate them itself.
# # So, first I need to create a VCF file. I am going to use freebayes. I want to direct freebayes to only the sites that have passed
# # the filters so far. It only takes intervals as input. So I first need to convert the keep.sites file to BED format.
# # I have written a Perl script called sites2bed.pl that turns a "sites" file, as created above for ANGSD, into a BED file, 
# # i. e. it combines consecutive sites into a BED interval. I have then used this script to turn the file of filtered sites
# # into a BED file, tha I can use for freebayes:
# ./sites2bed.pl ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites > ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.bed
# #
# # next I am splitting this BED file to allow parallisation of freebayes:
# mkdir split_beds
# split -l 1000 ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.bed split_beds/
# #
# # I am now running freebayes on the slimmed BAM files, parallelising over the split BED files:
# ls split_beds/* | 
# 	parallel -k -j 20 "freebayes -L slim.bamfile.list -f Big_Data_ref.fa -t {} --populations pop.list -T 0.01 -n 6 -m 2 -F 0.2 -V -j" \
# 	| vcffirstheader | vcfstreamsort | vcfuniq | bgzip > freebayes_ParEry.vcf.gz
# # 
# # Freebayes creates a version 4.1 VCF file, reporting GL (log10 scaled genotype likelihoods). 
# # I want to use 'vt' (http://genome.sph.umich.edu/wiki/Vt) to estimate the inbreeding coefficient F for each variable site based on
# # genotype likelihoods. However, vt requires phred-scaled genotype likelihoods. Fortuneately, bcftools comes with a plugin 'tag2tag' that converts
# # GL to PL.
# bcftools +tag2tag freebayes_ParEry.vcf.gz -- --gl-to-pl -r | vt compute_features - | bgzip > freebayes_ParEry.Vt_features.vcf.gz
# # This adds several new tags to the INFO column of the VCF, among which FIC (GL based inbreeding coefficent) and HWE_LPVAL (log p-value for deviation from
# # HWE, i. e. F = 0). The calculation of the Fis from genotype likelihoods is described here:
# # http://genome.sph.umich.edu/wiki/Genotype_Likelihood_based_Inbreeding_Coefficient
# # I think, the way it is described there, Fis has a upper limit of 1 and NO lower limit. Therefore, Fis below -1 are possible and they can be found in the 
# # VCF produced from vt.
# #
# # I want to filter this VCF for sites that have a significantly negative Fis. 'vcffilter' cannot filter for negative numbers. So I am using bcftools instead:
# bcftools view -H --include "INFO/FIC < 0 && INFO/HWE_LPVAL < -3" freebayes_ParEry.Vt_features.vcf.gz | cut -f 1 | sort -V | uniq > neg_Fis.contigs
# # This filters for negative F and a log(p-value) of less than -3 (e^-3 ~ 0.05). This finds 66 contigs with significantly negative F values, indicating
# # the mapping of paralogous sequences.
# #
# # Now, I am going to remove those contigs from the keep.sites file for ANGSD, creating a new sites file in the process:
# perl -ne'chomp; $H{$_}=1; END{$f=<*sites>; open(I, $f); while(<I>){@line=split; print if not exists $H{$line[0]}}}' neg_Fis.contigs \
# 	> ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites
# # This new sites file now contains no sites from contigs where the programme vt detected a sign. negative inbreeding coefficient.
# # The new keep.sites file containes 1,819,820 sites.
# #
# #
# # -----------------------------------------
# # TRIAL RUN THROUGH ANGSD/NGSTOOLS TUTORIAL
# # -----------------------------------------
# cd /data3/claudius/Big_Data/ANGSD
# # index new filtered sites file
# angsd sites index ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites
# # create region file containing a unique list of chromosomes:
# cut -f 1 ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites | uniq >ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.rf
# #
# # Quality Control
# #
# angsd -P 4 -b slim.bamfile.list -ref Big_Data_ref.fa -out Results/ParEry.qc -only_proper_pairs 0 \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites -minQ 0 -minMapQ 5 \
# 	-doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 1000 \
# 	-rf ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.rf
# # this took 23 min to complete. Without the regions file (*.rf) it is very slow.
# # This produces global and per-sample depth distribution files. Each column in the depth files
# # contains the number of sites in the depth category, from 0 till >1000, i. e. 1002 columns.
# #
# # The ngsTools tutorial (https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md) provides a couple of R scripts
# # to analyse and plot the output of ANGSD:
# Rscript Scripts/plotQC.R Results/ParEry.qc
# # This produces, among other, a pdf with three plots: quality score distribution, global depth distribution and per sample depth distribution for
# # the categories 0 - 9. There still seem to be sites with total depth above 1000. I can filter them out with ANGSD's '-setMaxDepth' filter.
# #
# #
# # reduced reference: only containing filtered contigs
# # further slimmed BAMs: using SLIM.BAMs.pl that uses the new rf file (noNegFis) 
# #
# # PCA
# #
# angsd -b SLIM.bamfile.list -ref reduced_ref.fa -out Results/ParEry -only_proper_pairs 0 \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites -minQ 20 \
# 	-minMapQ 5 -setMaxDepth 900 -minInd 18 -doCounts 1 -GL 1 -doMajorMinor 1 -doMaf 1 \
# 	-skipTriallelic 1 -doGeno 32 -doPost 1 -minMaf 0.01388 -baq 1
# 
# N_SITES=`zcat Results/ParEry.mafs.gz | tail -n +2 | wc -l`
# # angsd reports 85,773 variable sites int he *mafs file.
# 
# ls split_noNegFis.sorted.rf/* | parallel -j 20 "angsd -b SLIM.bamfile.list -ref reduced_ref.fa -out Results_split/ParEry.Geno8.{/} -only_proper_pairs 0 \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites -minQ 20 -minMapQ 5 -setMaxDepth 900 -minInd 18 -doCounts 1 -GL 1 -doMajorMinor 1 \
# 	-doMaf 1 -skipTriallelic 1 -doGeno 8 -doPost 1 -minMaf 0.01388 -baq 1 -rf {}"
# 
# #
# # POPULATION DIFFERENTIATION
# #
# # calculate sample allele frequency (SAF) likelihoods for PAR and ERY separately:
# for POP in ERY PAR; 
# do 
# 	echo $POP; 
# 	angsd -bam $POP.SLIM.bamfile.list \
# 	-ref reduced_ref.fa -anc reduced_ref.fa -out Results/$POP.unfolded \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites \
# 	-rf ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.rf \
# 	-only_proper_pairs 0 -baq 1 -minMapQ 5 -minQ 20 -minInd 10 -doCounts 1 \
# 	-setMaxDepth 450 -GL 1 -doSaf 1 -doMajorMinor 1 -skipTriallelic 1 -nThreads 2; 
# done
# # Note, that I am using the reference sequence to polarise the SAF's.
# #
# # estimate ML unfolded SFS for PAR and ERY population by EM:
# realSFS PAR.unfolded.saf.idx -P 20 2> /dev/null > PAR.unfolded.sfs
# realSFS ERY.unfolded.saf.idx -P 20 2> /dev/null > ERY.unfolded.sfs
# #
# # estimate ML unfolded 2D-SFS:
# realSFS -P 12 ERY.unfolded.saf.idx PAR.unfolded.saf.idx 2> /dev/null > ERY.PAR.unfolded.2dsfs
# #
# # calculating Fst, using 2D-SFS as prior:
# realSFS fst index ERY.unfolded.saf.idx PAR.unfolded.saf.idx -sfs ERY.PAR.unfolded.2dsfs -fstout ERY.PAR -P 12
# # The 2D-SFS is the prior, which is multiplied by the joint sample allele frequency likelihoods and thus turns them
# # into joint SAF posterior probabilities (empirical Bayes), e. g. the prob. of k minor alleles in pop1 and z minor alleles
# # in pop2. For each combination of number of minor alleles in each pop, an Fst is calculated and weighted by the joint
# # SAF posterior probabilities. realSFS returns the integral (sum) of this distribution, that is the posterior expectation.
# # That means, realSFS actually does calculate the posterior probability distribution of Fst's (see equation 16 in Fumagalli2013). 
# perl -e'print "Chrom\tPos\ta\taplusb\n"' > ERY.PAR.fst.tab
# realSFS fst print ERY.PAR.fst.idx 2> /dev/null >> ERY.PAR.fst.tab
# # The tab file can now be read into R.
# #
# realSFS fst stats ERY.PAR.fst.idx
# # this prints genome-wide Fst averaged over positions
# # The FST.Unweight seems to be the average of ratios of a and a+b (i. e. per position Fst estimates).
# # The Fst.weight seems to be the ratio of averages of a and a+b, i. e. avg(a)/avg(a+b).
# # In the output of `realSFS fst print ERY.PAR.fst.idx`, 'a' is in the third column and 'a+b' is in the fourth column.
# # This is confirmed by an issue discussion on the ANGSD repo: https://github.com/ANGSD/angsd/issues/16
# # Use ratio of averages estimator (Weighted)! See Bhatia2013
# # If supplying ""-whichFst 1", then realSFS calculates Hudson's Fst as shown in Bhatia2013, equation (10).
# #
# # NUCLEOTIDE DIVERSITY 
# #
# # compute allele frequency posterior probabilities using the FOLDED SFS as prior:
# for POP in ERY PAR; 
# do 
# 	echo $POP; 
# 	angsd -b $POP.SLIM.bamfile.list -ref reduced_ref.fa -anc reduced_ref.fa -fold 1 \
# 	-out Results/$POP -only_proper_pairs 0 -baq 1 -sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.sites \
# 	-minMapQ 5 -minQ 20 -minInd 10 -doCounts 1 -setMaxDepth 450 -gl 1 -doSaf 1 -domajorminor 1 -skipTriallelic 1 -doThetas 1 \
# 	-pest Results/$POP.sfs -rf ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.rf; \
# done 
# # This produces a *theta.gz file for each population. Since the output format is in textual, this command could also be parallelised over regions.
# #
# # Now let's calculate Watterson's theta, pi (Tajima's theta) and Tajima's D for each contig:
# # index thetas.gz files:
# thetaStat make_bed ERY.thetas.gz
# thetaStat make_bed PAR.thetas.gz
# # compute diversity statistics per contig:
# thetaStat do_stat ERY.thetas.gz -nChr 18
# thetaStat do_stat PAR.thetas.gz -nChr 18
# # Note, -nChr takes #ind when *thetas.gz was calculated from a folded SFS and 2*#ind when calculated from an unfolded SFS.
# # The output files *pestPG contain a few more statistics, but with a folded SFS only tW, tP and Tajima are sensible, the others cannot be interpreted without 
# # an unfolded SFS.
# 
# # # realSFS fst print EryPar.fst.idx 2> /dev/null | \
# # # 	perl -ne'BEGIN{print "chrom\tpos\ta\ta+b\ta/(a+b)\n"}{@l=split;chomp;print "$_\t";if($l[3]>0){print $l[2]/$l[3], "\n"}else{print "0\n"}}' > EryPar.fst.tab
# 
# 
# ################################## RESTART ANGSD ANALYSIS ##########################################
# #
# # HWE with ANGSD again
# #
# mkdir SnpStat
# #
# # estimate total sample F:
# #
# ls split_rf/* | \
# 	parallel -j 12 "angsd -b slim.bamfile.list -doSnpStat 1 -doMaf 1 -domajorminor 1 -skipTriallelic 1 -gl 1 -snp_pval 1e-2 \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SnpStat/{/} -rf {}"
# # Note, that I am _not_ using '-hwe_pval' here and I am filtering for sites with some evidence for being polymorphic via '-snp_pval'.
# # see ANGSD issue: https://github.com/ANGSD/angsd/issues/55
# #
# cd SnpStat
# # combine split output files:
# zcat aa.hwe.gz | head -n 1 > ParEry.hwe
# for f in *hwe.gz; do zcat $f | tail -n +2 >> ParEry.hwe; done
# # 
# # filter for sites with negative F and p-value below 0.05:
# gawk '$7<0 && $9 < 0.05' ParEry.hwe | cut -f 1 | uniq > ParEry.negFis.contigs
# # this finds 81 contigs
# # only 41 of these were also detected by vt from freebayes likelihoods (see above):
# grep -f ParEry.negFis.contigs ../neg_Fis.contigs | wc -l
# #
# # create new keep.sites file:
# cd ..
# perl -ne'chomp; $H{$_}=1; 
# 	END{open(I, "ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites"); 
# 	while(<I>){@line=split; print if not exists $H{$line[0]}}}' SnpStat/negFis.contigs \
# 	> ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSDnegFisFiltered.sorted.sites
# # This new keep.sites file contains 1,817,011 on 34,886 contigs.
# #
# # estimate within population F
# #
# for POP in ERY PAR; 
# do 
# 	ls split_rf/* | \
# 	parallel -j 12 "angsd -bam $POP.slim.bamfile.list -dosnpstat 1 -domaf 1 -domajorminor 1 -skiptriallelic 1 -gl 1 \
# 	-snp_pval 1e-2 -sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SnpStat/$POP.{/} -rf {}"; 
# done
# # This finds 180 contigs for ERY and only 41 for PAR. Only 15 of the contigs detected in PAR were also detected in ERY.
# # So, in total the population-wise HWE filter detected 206 contigs, compared to 81 contigs when using both populations together.
# # 70 contigs were detected by both the total sample HWE filter AND the population-wise HWE filter. That means that the population-wise
# # HWE filter detected an additional 136 contigs compared to the total sample HWE filter. 11 contigs were found only by the total sample
# # HWE filter. tHAT MEANS, BOTH hwe FILTERS TOGETHER HAVE DETECTED 217 CONTIGS:
# cat PAR.negFis.contigs ERY.negFis.contigs ParEry.negFis.contigs | sort -V | uniq >combined.negFis.contigs
# #
# # create new keep.sites file:
# #
# cd /data3/claudius/Big_Data/ANGSD
# perl -ne'chomp; $H{$_}=1;\
# 	END{open(I, "ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.sorted.sites");while(<I>){@line=split; print if not exists $H{$line[0]}}}' \
# 	SnpStat/combined.negFis.contigs \
# 	> ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.sites
# # THE NEW KEEP.SITES FILE CONTAINS 1,799,962 SITES ON 34,750 CONTIGS.
# # contigs with F close to one could be affected by allele dropout (due to polymorphism in the restriction site)
# # or they could map to the X chromosome, since all sequenced individuals are males 
# 
# #
# # Quality Control
# #
# cd /data3/claudius/Big_Data/ANGSD
# mkdir Quality_Control
# #
# # Note, the following ANGSD command is very slow. Use the samtools depth command right below it instead.
# #
# # Note, I cannot use reduced_ref.fa produced above for this command, because ANGSD throughs an error saying that it can't find
# # a contig in the reference that was in the BAM. That is because reduced_ref.fa was produced using keep.rf file produced with HWE filtering
# # by vt. Since the BAM files referred to in the slim.bamfile.list contain all reads (and a few more)
# # from the regions in the *rf file, this command should do the right thing.
# angsd -P 4 -b slim.bamfile.list -ref Big_Data_ref.fa -out Quality_Control/ParEry.qc -only_proper_pairs 0 \
# 	-sites ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.sites \
# 	-minQ 0 -minMapQ 5 -doDepth 1 -doCounts 1 -dumpcounts 2 \
# 	-rf ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.rf
# 	-maxDepth 1000
# # I have checked that *counts.gz file reports depths for all positions in the keep.sites file, i. e. 1,799,962.
# # There shouldn't be positions with total depth less than 3x15=45, since that is what I filtered for with 'even.depth.pl'.
# cd /data3/claudius/Big_Data/ANGSD
# # I used sites2bed.pl to convert the keep.sites file into bed format for the following command:
# samtools depth -aa -b ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.bed -Q 5 Data/*sorted.slim.bam | \
# 	bgzip > Quality_Control/ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.depths.gz
# # Note that there 12,693 positions with coverage greater than 1000x.
# #
# #
# # EXCESS ACROSS-SAMPLE COVERAGE
# #
# # Get the 99th percentile of the global (i. e. across sample) coverage distribution:
# cd /data3/claudius/Big_Data/ANGSD/Quality_Control
# zcat ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.depths.gz | \
# 	perl -ane'BEGIN{use List::Util qw(sum)} print "@F[0..1]", "\t",sum(@F[2..$#F]), "\n"' | \
# 	gawk '{print $3}' | \
# 	Rscript -e 'x=scan("stdin",quiet=T); cat(quantile(x, probs=c(0.99)),fill=T)'
# 	> global.Q99.cov
# # the global Q99 is 914
# #
# # the following command creates a list of contigs with positions that have total coverage above the global Q99:
# zcat ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.depths.gz | \
# 	perl -ane'BEGIN{use List::Util qw(sum)} print "@F[0..1]", "\t",sum(@F[2..$#F]), "\n"' | \
# 	gawk '$3>914' | cut -d" " -f 1 | uniq > gtGlobalQ99.contigs
# # This finds 407 additonal contigs that should be exluded from further analysis due to excessive global coverage.
# #
# # create new keep.sites file:
# #
# cd ..
# perl -ne'chomp; $H{$_}=1;
# 	END{open(I, "ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.sorted.sites"); 
# 	while(<I>){@line=split; print if not exists $H{$line[0]}}}' Quality_Control/gtGlobalQ99.contigs \
# 	> ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtGlobalQ99Cov.sorted.sites
# # This new keep.sites file contains 1,730,524 sites on 34,343 contigs.
# #
# # After converting the new keep.sites file into bed format with sites2bed.pl, I can run the samtools depth command again in order
# # to get the new coverage distributions:
# samtools depth -aa -b ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted.bed \
# 	-Q 5 Data/*sorted.slim.bam | \
# 	bgzip > Quality_Control/ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted.depths.gz 
# #
# #
# # ------------------------------------- PCA ---------------------------------------------------------------------------------------------
# #
# # With global MAF estimate as prior for genotype probabilities
# #
# # KNOWN MINOR ALLELE and SNP CALLING:
# #
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# mkdir SPLIT_RF
# split -l 500 $keep.rf SPLIT_RF/
# ls SPLIT_RF/* | parallel -j 12 "angsd -rf {} -bam slim.bamfile.list -ref Big_Data_ref.fa -out PCA/GlobalMAFprior/EryPar.{/} \
# 	-only_proper_pairs 0 -sites $keep.sites -minMapQ 5 -baq 1 -doCounts 1 -GL 1 -domajorminor 1 -doMaf 1 -skipTriallelic 1 \
# 		-SNP_pval 1e-3 -doGeno 32 -doPost 1"
# # Note, that the assumption of HWE is violated here, but that is how Matteo shows it for 3 divergent human populations (see ngsTools tutorial)
# # and Stefan does it for his Hawaiian spiders (see Ariamnes_code.pdf)
# # It also seems that without specifying -minInd, the *mafs file contains only sites where there is data from all 36 individual.
# # I have gunzip'ed all *geno.gz files, which are binary, and concatenated them simply with 'cat'.
# # Then I counted the number of sites for which MAF's are reported from the *mafs files, which are textual.
# # ANGSD reports MAF's and derived from this posterior genotype probabilities for 73,841 variable sites (with a p-value of being variable
# # less than 1e-3).
# ngsCovar -probfile EryPar.geno -nind 36 -nsites 73841 -call 0 -norm 0 -outfile ParEry.covar
# # I am not using called genotypes and I am not normalising by allele frequency here.
# #
# # ###
# # PCA with UNKNOWN MINOR ALLELE (-doMaf 2) and WITHOUT snp CALLING
# # ###
# # I want to redo the covariance matrix calculation but this time I want to take uncertainty in the determination of the minor allele frequency 
# # into account by using -doMaf 2 (instead -doMaf 1) and I want to replace the SNP calling by providing a file containing posterior SAF's
# # to ngsCovar (see Fumagalli2013).
# #
# # 1) calcualte genotype posterior probabilities using MAF's as prior (and HWE assumption):
# cd /data3/claudius/Big_Data/ANGSD
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# mkdir PCA/GlobalMAFprior/noSNPcall_unknownMinor
# ls SPLIT_RF/* | parallel -j 12 "angsd -rf {} -bam slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa \
# 	-out PCA/GlobalMAFprior/noSNPcall_unknownMinor/EryPar.{/} \
# 	-only_proper_pairs 0 -sites $keep.sites -minMapQ 5 -baq 1 -GL 1 -domajorminor 1 -doMaf 2 \
# 	-doGeno 32 -doPost 1"
# # Note, merging of *saf files with 'realSFS cat' doesn't work on huluvu (see https://github.com/ANGSD/angsd/issues/60).
# # So, I have to run SAF calculation without GNU parallel.
# #
# # 2) calculate SAF posterior probabilities
# #
# # Sample Allele frequency likelihoods UNFOLDED
# #
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# mkdir -p SAFs/ERY SAFs/PAR SAFs/EryPar
# # GLOBAL:
# angsd -rf $keep.rf -bam slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa \
# 	-out SAFs/EryPar/EryPar.unfolded \
# 	-only_proper_pairs 0 -sites $keep.sites -minMapQ 5 -baq 1 -GL 1 -doSaf 1 -fold 0 -P 1
# # Note, no minimum individual filtering, I want SAF's for every site. 
# # There are 1,730,389 sites in the *saf.idx file, 135 fewer than in the keep.sites file, but the same number as
# # in the mafs output from the upper angsd -domaf 2 command.
# # Turning on multithreading, makes the upper command very slow. Therefore, -P 1.
# #
# # 3) determine global unfolded SFS
# #
# # GLOBAL UNFOLDED SFS:
# cd /data3/claudius/Big_Data/ANGSD/SAFs/EryPar
# realSFS EryPar.unfolded.saf.idx > EryPar.unfolded.sfs
# #
# # 4) calculate posterior probabilities of SAF's
# #
# # Now, I am using this global unfolded SFS as prior for another ANGSD run in order to calculate posterior probabilities of SAF's:
# cd /data3/claudius/Big_Data/ANGSD
# angsd -rf $keep.rf -bam slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa \
# 	-out SAFs/EryPar/EryPar.unfolded.postProb \
# 	-only_proper_pairs 0 -sites $keep.sites -minMapQ 5 -baq 1 -GL 1 -doSaf 1 -fold 0 -P 1 \
# 	-pest /data3/claudius/Big_Data/ANGSD/SAFs/EryPar/EryPar.unfolded.sfs
# # I unzip'ed the *saf.gz file.
# # I created a symbolic link to this *saf file in /data3/claudius/Big_Data/ANGSD/PCA/GlobalMAFprior/noSNPcall_unknownMinor.
# #
# # 5) use SAF posterior probabilities to determine P_var,s (eq. 20 in Fumagalli2013) in order to weight each site by its
# #    probability of being variable during the calculation of covariances:
# cd /data3/claudius/Big_Data/ANGSD/PCA/GlobalMAFprior/noSNPcall_unknownMinor
# ngsCovar -probfile EryPar.geno -outfile EryPar.covar -nind 36 -nsites 1730389 -call 0 -norm 0 -sfsfile EryPar.unfolded.postProb.saf
# # ###
# #
# #
# # PCA with unknown minor allele (-doMaf 2) and SNP calling:
# #
# mkdir /data3/claudius/Big_Data/ANGSD/PCA/GlobalMAFprior/withSNPcall_unknownMinor
# cd /data3/claudius/Big_Data/ANGSD
# ls SPLIT_RF/* | parallel -j 12 "angsd -rf {} -bam slim.bamfile.list -ref Big_Data_ref.fa -out PCA/GlobalMAFprior/withSNPcall_unknownMinor/EryPar.{/} \
# 	-only_proper_pairs 0 -sites $keep.sites -minMapQ 5 -baq 1  -GL 1 -domajorminor 1 -doMaf 2 \
# 		-SNP_pval 1e-3 -doGeno 32 -doPost 1"
# # There are 68,590 sites reported in the *mafs files, i. e. that many were called as variant. 
# ngsCovar -probfile EryPar.geno -nsites 68590 -nind 36 -call 0 -norm 0 -outfile EryPar.covar
# #
# #
# # PCA with unknown minor allele (-doMaf 2) and SNP calling and genotype calling:
# #
# # I am using the same output as from the previous section and use ngsCovar to do the genotype calling
# # from the genotype posterior probabilities (in the *geno file). The genotype with the maximum posterior
# # probability is called:
# ngsCovar -probfile EryPar.geno -nsites 68590 -nind 36 -call 1 -norm 0 -outfile EryPar.covar.GC
# 
# 
# 
# # ------------------------------------- FST ---------------------------------------------------------------------------------------------
# 
# # 1) calculate UNFOLDED SAF's   
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# mkdir -p SAFs/EryPar
# # PER POP:
# angsd -bam PAR.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/PAR/PAR.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1
# # this calculates SAF's for 1,214,939 sites. I required at least 9 individuals to have reads. That's why it is less than 1.7M.
# angsd -bam ERY.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/ERY/ERY.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1
# # this calculates SAF's for 1,638,468 sites.
# 
# # 2) estimate ML unfolded 2D-SFS:
# mkdir FST
# realSFS -P 12 SAFs/ERY/ERY.unfolded.saf.idx SAFs/PAR/PAR.unfolded.saf.idx 2> /dev/null > FST/EryPar.unfolded.2dsfs
# # It is important to understand the output format of the *.2dsfs file produced by realSFS:
# # Here, I have first given the SAF's from Ery, then the SAF's from Par to realSFS. The *.2dsfs file contains just 
# # a single line of numbers, i. e. the flattened 2D matrix. The matrix has a ERY x PAR dimension of 37x37, since
# # there were 18 diploid individuals in each population sample. The first 37 numbers in EryPar.unfolded.2dsfs contain
# # 2dsfs[1, 1:37], i. e. the counts of non-reference alleles that have count of 0 in ERY and 0-36 in PAR. The next 37
# # numbers are 2dsfs[2, 1:37], i. e. the counts of non-reference alleles that have count of 1 in ERY and 0-36 in PAR,
# # and so on.
# #
# tr ' ' '\n' < FST/EryPar.unfolded.2dsfs | gawk '{sum+=$1}END{print sum}'
# # this reports 1.13 M counts in the 2D-SFS. So, the  ML unfolded 2D-SFS was calculated from 1.13 M sites,
# # that had reads from at least 9 individuals in each population. realSFS has automatically determined the overlapping sites
# # between the two populations.
# 
# # 3) use unfolded 2D-SFS as prior for Fst estimation
# # Reynold's Fst:
# realSFS fst index SAFs/ERY/ERY.unfolded.saf.idx SAFs/PAR/PAR.unfolded.saf.idx -sfs FST/EryPar.unfolded.2dsfs -fstout FST/EryPar.Reynolds -P 12
# # Hudson/Bhatia's Fst:
# realSFS fst index SAFs/ERY/ERY.unfolded.saf.idx SAFs/PAR/PAR.unfolded.saf.idx -sfs FST/EryPar.unfolded.2dsfs -whichFst 1 -fstout FST/EryPar.Bhatia -P 12
# #
# realSFS fst print EryPar.Bhatia.fst.idx 2> /dev/null > EryPar.Bhatia.fst.tab
# realSFS fst print EryPar.Reynolds.fst.idx 2> /dev/null > EryPar.Reynolds.fst.tab
# # These two files each contain per-site numerators (col 3) and denominators (col 4) of Hudson/Bhatia's or Reynold's Fst, respectively.
# # For Hudson/Bhatia's Fst see eq. (9) and (10) in Bhatia2013. For Reynold's Fst see eq. (1)-(3) in Fumagalli2013.
# # also see ANGSD issue https://github.com/ANGSD/angsd/issues/61
# # There are Fst's for 1,629,161 sites in the *fst.idx files, but only 1.13 M sites have reads from 9 or more individuals in both populations.
# # ==> create sites files with different minInd thresholds for each population, use this to filter sites from Fst output, calculate average Fst
# # for different minInd thresholds, see Gautier2013.
# # There 32,706 contigs in each *fst.idx file.
# 
# # ---- FOLDED SPECTRA ----
# #
# # see https://groups.google.com/d/msg/ngstools-user/9Z4viLJ7NJA/RWfQyoXnDgAJ
# #
# # 1) calculate FOLDED SAF's
# # PER POP:
# angsd -bam PAR.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/PAR/PAR.FOLDED -fold 1 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1
# # this calculates SAF's for 1,214,939 sites. I required at least 9 individuals to have reads. That's why it is less than 1.7M.
# angsd -bam ERY.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/ERY/ERY.FOLDED -fold 1 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1
# # this calculates SAF's for 1,638,468 sites.
# #
# # 2) estimate ML FOLDED 2D-SFS:
# realSFS -P 12 SAFs/ERY/ERY.FOLDED.saf.idx SAFs/PAR/PAR.FOLDED.saf.idx 2> /dev/null > FST/EryPar.FOLDED.2dsfs
# 
# 
# # 3) use FOLDED 2D-SFS as prior for Fst estimation
# # Reynold's Fst:
# realSFS fst index SAFs/ERY/ERY.FOLDED.saf.idx SAFs/PAR/PAR.FOLDED.saf.idx -sfs FST/EryPar.FOLDED.2dsfs -whichFst 0 -fstout FST/EryPar.FOLDED.Reynolds -P 12
# # Hudson/Bhatia's Fst:
# realSFS fst index SAFs/ERY/ERY.FOLDED.saf.idx SAFs/PAR/PAR.FOLDED.saf.idx -sfs FST/EryPar.FOLDED.2dsfs -whichFst 1 -fstout FST/EryPar.FOLDED.Bhatia -P 12
# #
# realSFS fst print EryPar.FOLDED.Bhatia.fst.idx 2> /dev/null > EryPar.FOLDED.Bhatia.fst.tab
# realSFS fst print EryPar.FOLDED.Reynolds.fst.idx 2> /dev/null > EryPar.FOLDED.Reynolds.fst.tab
# 
# # --- BOOTSTRAP ----
# #
# # I want to randomise the population labels of individuals and calculate Fst for each randomisation
# #
# # 1) create randomised bamfile.list files
# for i in {1..100}; do shuf slim.bamfile.list | tee >(head -n 18 >SAFs/Bootstrap/Pop1/$i.bamfile.list) | tail -n 18 > SAFs/Bootstrap/Pop2/$i.bamfile.list; done
# #
# # 2) calculate unfolded SAF's for each randomisation
# cd /data3/claudius/Big_Data/ANGSD
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# # Pop1:
# seq 1 100  | parallel -j 12 "angsd -bam SAFs/Bootstrap/Pop1/{}.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/Bootstrap/Pop1/{}.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1"
# # Pop2:
# seq 1 100  | parallel -j 12 "angsd -bam SAFs/Bootstrap/Pop2/{}.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/Bootstrap/Pop2/{}.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1"
# #
# # 3) estimate 2D-SFS
# seq 1 100 | parallel -j 6 "realSFS -P 4 SAFs/Bootstrap/Pop1/{}.unfolded.saf.idx SAFs/Bootstrap/Pop2/{}.unfolded.saf.idx 2> /dev/null > FST/Bootstrap/{}.unfolded.2dsfs"
# # this takes a couple of hours
# #
# # 4) estimate Fst
# seq 1 100 | parallel -j 6 "realSFS fst index SAFs/Bootstrap/Pop1/{}.unfolded.saf.idx SAFs/Bootstrap/Pop2/{}.unfolded.saf.idx \
# 	-sfs FST/Bootstrap/{}.unfolded.2dsfs -whichFst 1 -fstout FST/Bootstrap/{}.unfolded.Bhatia -P 4 2>/dev/null"
# #
# # 5) extract Fst table
# cd FST/Bootstrap
# for f in *fst.idx; do realSFS fst print $f 2>/dev/null > `basename $f idx`tab; done
# #
# #
# #
# # ---- Fst for different ascertainment classes ----
# #
# # As in Bhatia2013, fig. 1, I want to classify sites by minor allele frequency in either PAR or ERY.
# # I do not need to recalculate Fst, I just need to estimate MAF's for each population separately.
# cd /data3/claudius/Big_Data/ANGSD
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# mkdir FST/MAFs/PAR FST/MAFs/ERY
# for POP in PAR ERY;
# do
# 	echo $POP;
# 	ls SPLIT_RF/* | parallel -j 12 "angsd -rf {} -bam $POP.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa \
# 	-out FST/MAFs/$POP/{/} \
# 	-only_proper_pairs 0 -sites $keep.sites -minMapQ 5 -baq 1 -GL 1 -domajorminor 1 -doMaf 2";
# 	sleep 10;
# done
# # There are 1,730,524 sites in the keep.sites file. 
# # For 1,728,745 of them, MAF's have been calculated in ERY
# # and for 1,710,005 of them MAF's have been calculated in PAR. 
# #
# # I need to extract from the MAF files those sites for which I have Fst's calculated. 
# cd /data3/claudius/Big_Data/ANGSD/FST/MAFs/ERY
# perl -ane'$H{$F[0]}{$F[1]}=[@F];
# 	END{open(I, "tail -n +2 ERY.mafs | "); 
# 	while(<I>){
# 		@line=split; if(exists $H{$line[0]}{$line[1]}){chomp; print; print "\t", join("	", @{$H{$line[0]}{$line[1]}}), "\n"}; 
# 		}
# 	}' /data3/claudius/Big_Data/ANGSD/FST/EryPar.Bhatia.fst.tab \
# 	> ERY.mafs.withFST
# # There are 1,513,856 sites with ERY mafs and Fst.
# #
# # # alternative way to get the overlap:
# # join -t" " <(sed 's/^\(\w*\)\t/\1:/' ../../EryPar.Bhatia.fst.tab | sort -bk1,1) <(sed 's/^\(\w*\)\t/\1:/' ERY.mafs | sort -bk1,1) | sed 's/:/\t/' > ERY.MAFS.withFST
# # there are issues with the correct way of sorting for the join command
# #
# cd /data3/claudius/Big_Data/ANGSD/FST/MAFs/PAR
# perl -ane'$H{$F[0]}{$F[1]}=[@F];
# 	END{open(I, "tail -n +2 PAR.mafs | "); 
# 	while(<I>){
# 		@line=split; if(exists $H{$line[0]}{$line[1]}){chomp; print; print "\t", join("	", @{$H{$line[0]}{$line[1]}}), "\n"}; 
# 		}
# 	}' /data3/claudius/Big_Data/ANGSD/FST/EryPar.Bhatia.fst.tab \
# 	> PAR.mafs.withFST
# # There are 1,496,271 sites with PAR mafs and Fst.
# #
# # Both *withFST files can be read into R for analysis.
# #
# #
# #
# # ------------------------------ Nucleotide Diversity -----------------------------------
# #
# # I have previously calculated folded SAF's for Fst estimation (but that didn't work correctly with folded SAFs).
# # I am repeating the command lines that I used to estimate *folded* SAF likelihoods:
# # 1) calculate FOLDED SAF's
# # PER POP:
# angsd -bam PAR.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/PAR/PAR.FOLDED -fold 1 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1
# # this calculates SAF's for 1,214,939 sites. I required at least 9 individuals to have reads. That's why it is less than 1.7M.
# angsd -bam ERY.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/ERY/ERY.FOLDED -fold 1 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1
# # this calculates SAF's for 1,638,468 sites.
# #
# # I now want to estimate per population SFS's with bootstrapping.
# #
# # 2) estimate ML FOLDED SFS:
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 SAFs/ERY/ERY.FOLDED.saf.idx 2>/dev/null > SFS/ERY/ERY.FOLDED.sfs
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 SAFs/PAR/PAR.FOLDED.saf.idx 2>/dev/null > SFS/PAR/PAR.FOLDED.sfs
# # the upper two commands make sure that a good ML estimate is reached, albeit they reach it slowlier than the accelerated EM algorithm.
# # I've timed these commands: see m0.ery.time and m0.par.time
# 
# #
# # testing variation in ML estimate of SFS:
# #
# for i in {1..100};
# do
# 	realSFS -P 12 SAFs/ERY/ERY.FOLDED.saf.idx 2> /dev/null >> SFS/ERY/ERY.FOLDED.sfs.ml
# 	realSFS -P 12 SAFs/PAR/PAR.FOLDED.saf.idx 2> /dev/null >> SFS/PAR/PAR.FOLDED.sfs.ml
# done
# 
# #
# # 3) bootstrap ML FOLDED SFS:
# #
# realSFS -P 10 -bootstrap 1000 SAFs/ERY/ERY.FOLDED.saf.idx 2> /dev/null > SFS/ERY/ERY.FOLDED.sfs.boot
# realSFS -P 10 -bootstrap 1000 SAFs/PAR/PAR.FOLDED.saf.idx 2> /dev/null > SFS/PAR/PAR.FOLDED.sfs.boot
# #
# #
# # I am trying to do bootstrap resampling of SFS with exhaustive ML search:
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 -bootstrap 200 SAFs/ERY/ERY.FOLDED.saf.idx 2> /dev/null > SFS/ERY/ERY.FOLDED.sfs.boot.exh
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 -bootstrap 200 SAFs/PAR/PAR.FOLDED.saf.idx 2> /dev/null > SFS/PAR/PAR.FOLDED.sfs.boot.exh
# # this takes a while to finish
# 
# #
# # thetaStat
# #
# cd /data3/claudius/Big_Data/ANGSD
# mkdir -p THETASTAT/ERY THETASTAT/PAR
# export keep="ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted"
# #
# # calculate posterior folded SAF's and estimate posterior diversity statistics (-doThetas) with folded SFS as prior
# #
# # ERY
# ls SPLIT_RF/* | parallel -j 12 "angsd -b ERY.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -fold 1 \
# -out THETASTAT/ERY/ERY.{/} -only_proper_pairs 0 -baq 1 -sites $keep.sites -minMapQ 5 -minInd 9 -gl 1 -doSaf 1 -doThetas 1 \
# -pest SFS/ERY/ERY.FOLDED.sfs -rf {}"
# #
# # PAR
# ls SPLIT_RF/* | parallel -j 12 "angsd -b PAR.slim.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -fold 1 \
# -out THETASTAT/PAR/PAR.{/} -only_proper_pairs 0 -baq 1 -sites $keep.sites -minMapQ 5 -minInd 9 -gl 1 -doSaf 1 -doThetas 1 \
# -pest SFS/PAR/PAR.FOLDED.sfs -rf {}"
# # 
# # I concatenated the *thetas.gz files into an ERY.thetas.gz and PAR.thetas.gz file.
# #
# # Now let's calculate Watterson's theta, pi (Tajima's theta) and Tajima's D for each contig:
# # index thetas.gz files:
# thetaStat make_bed ERY.thetas.gz
# thetaStat make_bed PAR.thetas.gz
# # Output in the thetas.gz are the log scaled per site estimates of the thetas
# # compute diversity statistics per contig:
# thetaStat do_stat ERY.thetas.gz -nChr 18
# thetaStat do_stat PAR.thetas.gz -nChr 18
# # Note, -nChr takes #ind when *thetas.gz was calculated from a folded SFS and 2*#ind when calculated from an unfolded SFS.
# # The output files *pestPG contain a few more statistics, but with a folded SFS only tW, tP and Tajima are sensible, 
# # the others cannot be interpreted with a folded SFS.
# # thetaStat do_stat produces a *pestPG file with diversity estimators per contig.
# # Output in the pestPG file are the sum of the per site estimates for a contig.
# 
# # ---------------------------------------
# # Does folding need to be done by ANGSD ? 
# # ---------------------------------------
# # I want to compare folded SFS's with unfolded SFS's that have only been folded in dadi. Do they differ significantly?
# # I already created both folded and unfolded SAF's for each population:
# # + unfolded SAF's: lines 1423 - 1432
# # + folded SAF's: lines 1573 - 1582
# # I have already estimated folded SFS's for each population with realSFS: lines 1586 - 1590.
# # I still have to estimate unfolded 1D-SFS's for each population:
# 
# # 2) estimate ML unfolded 1D SFS:
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 SAFs/ERY/ERY.unfolded.saf.idx 2>/dev/null > SFS/ERY/ERY.unfolded.sfs
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 SAFs/PAR/PAR.unfolded.saf.idx 2>/dev/null > SFS/PAR/PAR.unfolded.sfs
# # the upper two commands make sure that a good ML estimate is reached, albeit they reach it slowlier than the accelerated EM algorithm.
# 
# 
# 
# 
# 
# 
# 
# ##############################################################################################
# # remove PCR duplicates
# ##############################################################################################
# # 28/04/2017
# 
# cd /data3/claudius/Big_Data
# 
# # note, it is necessary to first decrypt the input fastq files in the data directory
# 
# OUT=/data3/claudius/Big_Data/deduplicated
# 
# for se in data/*fq_1.gz;
# do
# 	pe=$(echo $se | sed 's/fq_1/fq_2/')
# 	out_pe=$(echo $pe | sed 's/\(^.*\)\(\.fq_2.*\)/\1_dedup\2/')
# 	out_se=$(echo $se | sed 's/\(^.*\)\(\.fq_1.*\)/\1_dedup\2/')
# 	starcode -d 2 -t 12 --connected-comp --non-redundant -1 <(zcat $se) -2 <(zcat $pe) --output1 >(gzip > $OUT/`basename $out_se`) --output2 >(gzip > $OUT/`basename $out_pe`)
# done
# # The --non-redundant option to starcode is different from the normal clustering ooutput format of starcode in that it keeps
# # some information about the input sequences, i. e. when given fastq input files it outputs fastq output files, instead of a
# # file containing just sequence strings. The starcode command collapses all read pairs that have a total edit distance across
# # single-end and paired-end reads of up to 2.
# 
# 
# # -----------------------------
# # count collapsed sequences
# # -----------------------------
# 
# cd /data3/claudius/Big_Data
# 
# sum=0; 
# 
# for file in deduplicated/*fq_1.gz;
# do 
# 	echo -n "$file   " >> deduplicated/uniqseq.count;
# 	count=`zcat $file | gawk '(NR-1)%4==0' | wc -l`;
#        	sum=$[$sum + $count];
#        	echo $count >> deduplicated/uniqseq.count;
# done 
# 
# echo "sum   $sum" >> deduplicated/uniqseq.count
#
# # The total number of retained read pairs is 10,661,569. This is very close to the number of fragments I previously counted from the same data set with my
# # Perl script purge_PCR_duplicates.pl.
# # Note, that starcode has a bug that turns on message-passing clustering in --non-redundant mode even when --connected-component clustering is specified.
# # That is the reason why the number of fastq sequences in the output is not identical to the output of my previous starcode clustering of read pairs
# # for de novo assembly. See the issue on github: https://github.com/gui11aume/starcode/issues/16



# # ----------------------------------------------------------
# # map deduplicated reads against Big Data reference 
# # ----------------------------------------------------------
# cd /data3/claudius/Big_Data/BOWTIE2
# mkdir BAM_dedup
# for SE in ../deduplicated/*fq_1.gz; do
# 	PE=`echo $SE | sed 's/fq_1/fq_2/'`
# 	ID=`basename $SE .fq_1.gz` 
# 	bowtie2 -x Big_Data_ref \
# 		-1 <(zcat $SE) -2 <(zcat $PE) \
# 		-q \
# 		--phred64 \
# 		--very-sensitive \
# 		--dpad 10 \
# 		--gbar 4 \
# 		--norc \
# 		--end-to-end \
# 		--np 10 \
# 		-I 60 \
# 		-X 800 \
# 		--no-contain \
# 		-t \
# 		--no-unal \
# 		--rg-id $ID \
# 		--rg "SM:$ID" \
# 		--omit-sec-seq \
# 		-p 22 \
# 		| samtools view -bq 1 - \
# 		> BAM_dedup/$ID.bam;
# done
# cd ..
# # this takes only about 15 min
# # running bowtie2 in end-to-end mode (i. e. no clipping of query sequences), "very sensitive" specifies a "seed" (kmer) length of 20
# # that needs to match exactly, seeds are sampled from the query according to the function f(x) = 1 + 0.5 * sqrt(x), x is the read length
# # with "--dpad 10" I am allowing gaps to be up to 10 bp long; "--gbar 4" disallows gaps within 4 bp of either end of the read; with
# # "--norc" bowtie2 will not attempt to align unpaired reads against the revcomp reference strand - I generally don't expect any good alignments
# # against the revcomp of the RADome, so this flag somewhat reduces unnecessary searches; "--np 10" sets the penalty for ambiguous 
# # characters in an alignment - in the RADome SE RAD tags and PE contigs are separated by up to 10 N's, I don't want alignments across this gap;
# # I specify a minimum fragment length of 60 (i. e. allowing some overlap between SE and PE reads) and a maximum fragment length of 800 (longer
# # than the longest sequence in the RADome; I am supressing the output of SAM records for reads that failed to align anywhere; if secondary sequences
# # are printed (which they shouldn't with this configuration), then omit the seq and qual string; "-p 22" runs bowtie2 on 22 cores; samtools view
# # then filters out all SAM records with 0 mapping quality and writes out a compressed BAM file


# # -------------------------------------------------
# # sort BAM files according to position in parallel
# # -------------------------------------------------
# parallel 'samtools sort -o {.}.sorted.bam {}' ::: ery*bam
# parallel 'samtools sort -o {.}.sorted.bam {}' ::: par*bam

# # -------------------------------------
# # change RG tag: add SM, LB and PL tag
# # -------------------------------------
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# rename -v 's/bam/bamx/' *sorted.bam
# for f in *bamx; do 
# 	n=`echo $f | sed 's/x$//'`;
# 	samtools view -h $f | \
# 	sed -r 's/^@RG\tID:(.*?)\t.*/@RG\tID:\1\tSM:\1\tLB:standard_RAD\tPL:GAIIx/' | \
#        	samtools view -b > $n; 
# 	rm -f $f;
# done 


# # -----------------------
# # mismappings of SE reads
# # -----------------------
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# samtools cat -o >(samtools view -f64 - | gawk '$4>2' | cut -f 3 | sort | uniq > exclude_contigs.SEposgt2) *bam
# # With this command line I am creating a list of contig ID's where a SE read from any ind mapped to a position greater than 2 
# # (remember, I padded the contigs with 1 N at beginning and end). A mapping position > 2 is a mismapping and indicates an
# # incorrect assembly of the contig. The whole contig should therefore be taken out of any further analysis. 
# # exclude_contigs.SEposgt2 contains the id's of 3632 contigs.


# # create a BED file noSEgt2.noDUST.bed
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# I have slighly modified the previous subtract.pl file.
# ./subtract_1.pl | sort -Vk1,1  > Big_Data_Contigs.noSEgt2.bed
# bedtools subtract -a Big_Data_Contigs.noSEgt2.bed -b ../Big_Data_ref.dust_intervals.sorted.bed > Big_Data_Contigs.noSEgt2.noDUST.bed
# 
# # exclude remainder of restriction site
# # From the 2nd till 7th position of each contig, positions cannot be polymorphic (TGCAGG). Inclusion of those sites would distort diversity estimates.
# # So I am going to remove those positions from the bed intervals:
# perl -lane'if($F[1] == 1){$F[1]=7; print join("  ", @F)}else{print}' Big_Data_Contigs.noSEgt2.noDUST.bed | gawk '$2 < $3' > Big_Data_Contigs.noSEgt2.noDUST.noTGCAGG.bed
# # note, the gawk command is necessary to make sure that no intervals are created where the starting position is greater or equal to the end position
  
# # create depth file with `samtools depth` 
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# samtools depth -aa -b Big_Data_Contigs.noSEgt2.noDUST.noTGCAGG.bed -Q 5 *sorted.bam | bgzip > ParEry.noSEgt2.noDUST.noTGCAGG.depth.gz
# # # This will print out a line for each position within the BED intervals. First column is contig name, second is position, the remaining
# # # 36 columns contain individual coverage for the site, from reads with mapQ > 5.
# 
# # # check that the number of sites in the .depth.gz file is the same as the length of the intervals in the BED file.
# gawk '{sum+=($3-$2)}END{print sum}' Big_Data_Contigs.noSEgt2.noDUST.noTGCAGG.bed
# zcat ParEry.noSEgt2.noDUST.noTGCAGG.depth.gz | wc -l
# # both commands return 85,488,084

# # ----------------------------
# # coverage filter
# # ----------------------------
# # I have written a python programme called coverage_filter.py (tested in coverage_filter.ipynb) that takes the samtools depth table
# # and filters out sites (not whole contigs) with excess coverage and not enough coverage. It defaults to filtering out sites with:
# #  - global coverage greater than the 99th percentile of the global coverage distribution
# #  - individual coverage greater than the 99th percentile of the individual coverage distribution (for each individual)
# #  - less than 15 individuals with at least 1x coverage
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# 
# /usr/bin/time -vo coverage_filter.time ./coverage_filter.py ParEry.noSEgt2.noDUST.noTGCAGG.depth.gz | gzip > ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.gz
# # this command took 3h 30min to finish
# # This filtering retained 768,439 sites from 85,488,084 sites (0.9%).
# zcat ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.gz | cut -f 1 | sort | uniq | wc -l
# # this prints 25,291 , i. e. sites from 25,291 contigs were retained.
# # Note, the previous analysis without removing PCR duplicates had 1,829,329 sites on 34,967 contigs at this stage of filtering. I have used a different
# # excess coverage filtering here: instead of discarding a whole contig if the SE read coverage was above the individual Q99, I am here discarding individual
# # sites if their individual coverage is above the respective individual Q99 of any individual or the global coverage is above the global Q99. I would have thought that this should
# # retain more sites not less. Previously, I had filtered for sites with at least 3x coverage in at least 15 individuals. Removing PCR duplicates and filtering
# # for at least 1x coverage in at least 15 individuals should not be a more stringent filter. So, it is not clear to me why removing PCR duplicates and coverage
# # filtering for sites not contigs did lead to less than half the number of sites retained by the previous analysis that used all reads, i. e. including PCR duplicates.

# 
# # --------------------
# # HWE filtering
# # --------------------
# # turn coverage filtered depth table into a sites file for ANGSD:
# zcat ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.gz | cut -f1,2 | sort -Vk1,1 -k2,2 > ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.sites
# 
# # continue analysis in the following directory
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# 
# # index sites file for ANGSD
# angsd sites index ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.sites
# 
# # create regions file for ANGSD
# cut -f 1 ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.sites | uniq > ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.rf
# 
# # split the regions file into small files, each containing only 500 contig names:
# split -l 500 ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.rf split_rf/

# # I am going to filter, as before, for negative deviations of genotype frequencies from HWE.
# # Large negative deviations should only be caused by mapping of paralogous sequences to the same position
# # in the reference.
# 
# # preparation:
# mkdir SnpStat
# mkdir data
# cd data
# ln -s ../../../BOWTIE2/BAM_dedup/*sorted.bam
# parallel samtools index {} ::: *bam
# cd ..
# 
# #
# # estimate total sample F:
# #
# ls split_rf/* | \
# 	parallel -j 12 "angsd -b bamfile.list -doSnpStat 1 -doMaf 1 -domajorminor 1 -skipTriallelic 1 -gl 1 -snp_pval 1e-2 \
# 	-sites ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SnpStat/{/} -rf {}"
# 
# # combine split output files:
# cd SnpStat
# zcat aa.hwe.gz | head -n 1 > ParEry.hwe
# for f in *hwe.gz; do zcat $f | tail -n +2 >> ParEry.hwe; done
# 
# #
# # filter for sites with negative F and p-value below 0.05:
# #
# gawk '$7<0 && $9 < 0.05' ParEry.hwe | cut -f 1 | uniq > ParEry.negFis.contigs
# # this finds 26 contigs
# 
# head -n 18 bamfile.list > ERY.bamfile.list
# tail -n 18 bamfile.list > PAR.bamfile.list
# 
# #
# # estimate within population F
# #
# for POP in PAR ERY; 
# do 
# 	ls split_rf/* | \
# 		parallel -j 12 "angsd -bam $POP.bamfile.list  -dosnpstat 1 -domaf 1 -domajorminor 1 -skiptriallelic 1 -gl 1 -snp_pval 1e-2 \
# 		-sites ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SnpStat/$POP.{/} -rf {}"; 
# done
# 
# # combine split output files:
# cd SnpStat
# for POP in ERY PAR;
# do
# 	zcat $POP.aa.hwe.gz | head -n 1 > $POP.hwe
# 	for f in $POP*hwe.gz; do zcat $f | tail -n +2 >> $POP.hwe; done
# done
# 
# #
# # filter out sites with negative F and p-value below 0.05:
# #
# for POP in ERY PAR;
# do
# 	gawk '$7<0 && $9 < 0.05' $POP.hwe | cut -f 1 | uniq > $POP.negFis.contigs;
# done
# # This finds 2 contigs in PAR and 35 in ERY.
# 
# # combine negFis.contigs files:
# cat PAR.negFis.contigs ERY.negFis.contigs ParEry.negFis.contigs | sort -V | uniq >combined.negFis.contigs
# # There are 50 contigs in the combined list of contigs
# 
# #
# # create new keep.sites file:
# #
# cd ..
# perl -ne'chomp; $H{$_}=1; \
# 	END{open(I, "ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.sites"); while(<I>){@line=split; print if not exists $H{$line[0]}}}' \
# 	SnpStat/combined.negFis.contigs > ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negFisFiltered.sites
# # I have checked that *negFisFiltered.sites contains 50 contigs fewer than *covFiltered.sites.
# 
# # Contigs with F close to one could be affected by allele dropout (due to polymorphism in the restriction site)
# # or they could map to the X chromosome, since all sequenced individuals are males. I am therefore going to use
# # WITHIN population Fis values to filter out contigs with an excess of homozygosity:
# #
# # filter positive Fis contigs
# #
# for POP in ERY PAR;
# do
# 	gawk '$7 > 0 && $9 < 0.05' $POP.hwe | cut -f 1 | uniq > $POP.posFis.contigs;
# done
# # this finds 220 contigs in ERY and 5 in PAR.
# 
# #
# # create new keep.sites file:
# #
# cd ..
# perl -ne'chomp; $H{$_}=1; \
# 	END{open(I, "ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negFisFiltered.sites"); while(<I>){@line=split; print if not exists $H{$line[0]}}}' \
# 	SnpStat/combined.posFis.contigs > ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negPosFisFiltered.sites
# # I have checked that the new keep.sites file contains 225 fewer contigs than the one before.
# # The new keep.sites file contains 756,502 sites on 25,016 contigs.

# # index sites file for ANGSD
# angsd sites index ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negPosFisFiltered.sites

# # create regions file for ANGSD
# cut -f 1 ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negPosFisFiltered.sites | uniq > ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negPosFisFiltered.rf
# 
# # -------------------------------------------------------
# # Sample Allele frequency likelihoods (SAF's) UNFOLDED
# # -------------------------------------------------------
# # # preparation:
# # mkdir -p SAFs/ERY SAFs/PAR 
# # ln -s ../Big_Data_ref.fa ../Big_Data_ref.fa.fai .
# 
# # PER POP:
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# export keep="ParEry.noSEgt2.noDUST.noTGCAGG.covFiltered.negPosFisFiltered"
# angsd -bam PAR.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/PAR/PAR.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 0 -GL 1 -doSaf 1 -nThreads 1 &
# angsd -bam ERY.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out SAFs/ERY/ERY.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 0 -GL 1 -doSaf 1 -nThreads 1 &
# # Note, no minimum individual filtering, I want SAF's for every site. 
# # running both commands in parallel (&) takes less than 50 min
# 
# cd SAFs/PAR
# realSFS print PAR.unfolded.saf.idx 2>/dev/null | less -S
# # The SAF files contain the loglikelihood ratio to the most likely allele frequency category.
# # https://github.com/ANGSD/angsd/blob/newsaf/doc/formats.pdf
# # The SAF file for PAR contains 748,141 sites on 24,770 contigs.
# # The SAF file for ERY contains 755,674 sites on 24,999 contigs.

# # --------------------------
# # estimate ML UNFOLDED SFS
# # --------------------------
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# mkdir -p SFS/PAR SFS/ERY
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 SAFs/ERY/ERY.unfolded.saf.idx 2>/dev/null > SFS/ERY/ERY.unfolded.sfs
# realSFS -P 12 -maxIter 50000 -tole 1e-6 -m 0 SAFs/PAR/PAR.unfolded.saf.idx 2>/dev/null > SFS/PAR/PAR.unfolded.sfs
# # this takes a about 2 hours to finish
# # the upper two commands make sure that a good ML estimate is reached, albeit they reach it slowlier than the accelerated EM algorithm.
# # === The resulting 1D spectra look terrible. They are completely unusable! ===
# # It could be that there simply is not enough information in the deduplicated read data to distinguish alleles from sequencing errors.
# # It is still not clear to me why my new coverage filtering with the deduplicated reads retained fewer sites on fewer contigs than my
# # previous coverage filtering with read data that included PCR duplicates and filtered whole contigs. I am therefore going to check
# # what influence the difference in coverage filtering had on the site frequency spectra.
# #
# # redo EXCESS COVERAGE FILTERING
# #
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# # In the following I am calculating SE coverage counts for each individual and
# # streaming those into R to determine the 99th percentile of the coverage distribution for each individual:
# echo "# filename       Q99_cov" > Q99.cov
# for f in *sorted.bam; 
# do 
# 	echo -n "$f " >> Q99.cov; \
# 	samtools view -f 64 $f | cut -f 3 | uniq -c | gawk '{print $1}' | \
# 	Rscript -e 'x=scan("stdin",quiet=T);cat(quantile(x,probs=c(0.99)),fill=T);' \
# 	>> Q99.cov;
# done
# # Note, that this inlcudes SE reads and their contigs that did not pass the mismapping filter above. 
# # In contrast to the depths in the output of the samtools depth command run on all *bam files (with -aa option),
# # this gets contigs (not sites) where the individual had at least one SE read mapped. Since the majority of sites
# # in the depth file from samtools have only reads mapped from one individual, the depth distribution from the samtools
# # depth command should have a mode of 0, while the distribution created by the above command line has a mode of 1 or above
# # (the 0 category doesn't even exist!).
# 
# # Next I want to use these 99 percentiles of the individual coverage distributions in order to detect contigs that have SE read coverage
# # above the individual Q99 in any individual. The following script does that.
# print_Q99_exCov_contigs.pl
# # This takes Q99.cov as input and spits out a list of contigs to exclude due to excessive coverage:
# # Big_Data_Contigs.gtQ99Cov
# # A contig has excess coverage if in any of the individual BAM files it has coverage above the 99th percentile
# # of that individual's coverage distribution. 
# # This finds 2192 contigs with excessive per individual coverage.
# 
# 
# # With the following command I create a list of all contigs that get at least one read mapped by any of the 36 individuals:
# samtools cat *sorted.bam | samtools view | cut -f 3 | uniq | sort | uniq > Big_Data_Contigs_with_mappings.list
# # The list contains 582,089 contigs (when using reads including PCR duplicates it were 582,596). 
# 
# #So, I have three files with contig names
# # # # to keep or delete from Big_Data_ref.fa.bed:
# # # # Big_Data_Contigs_with_mappings.list
# # # # Big_Data_Contigs.SEposgt2
# # # # Big_Data_Contigs.gtQ99Cov
# 
# # The Perl script subtract.pl takes these four files and outputs a new BED file containing intervals that passed the mismapping and excess
# # coverage filters:
# ./subtract.pl | sort -Vk1,1  > Big_Data_Contigs.noSEgt2.nogtQ99Cov.bed
# 
# # 583312 Big_Data_ref.fa.bed
# # 576463 Big_Data_Contigs.noSEgt2.nogtQ99Cov.bed
# #   6849 filtered out
# 
# # # # Now I need to subtract the DUST intervals from this BED file:
# bedtools subtract -a Big_Data_Contigs.noSEgt2.nogtQ99Cov.bed -b ../Big_Data_ref.dust_intervals.sorted.bed > Big_Data_Contigs.noSEgt2.nogtQ99Cov.noDUST.bed
# 
# # The total length of the new filtered intervals is:
# gawk '{sum+=($3-$2)}END{print sum}' Big_Data_Contigs.noSEgt2.nogtQ99Cov.noDUST.bed
# # 88,272,073
# 
# #
# # EVEN COVERAGE                                                                                                                            
# #
# # I want to filter individual sites for sufficiently even coverage across individuals. 
# samtools depth -aa -b Big_Data_Contigs.noSEgt2.nogtQ99Cov.noDUST.bed -Q 5 *sorted.bam | bgzip > ParEry.noSEgt2.nogtQ99Cov.noDUST.depth.gz
# 
# # Finally I am running a filter over this depth table, only keeping sites with 1 or more coverage in at least 15 of the individuals:
# ./even.depth.pl 1 15 ParEry.noSEgt2.nogtQ99Cov.noDUST.depth.gz
# # this only takes about 30 min
# # This should produce the file called:
# # ParEry.noSEgt2.nogtQ99Cov.noDUST.1.15.depth.gz
# # This file just contains two columns: contig name <TAB>  position
# # The positions are 1-based. There are 3,466,140 sites in this file (see line 1002 for comparison).
# 
# 
# # # From the 2nd till 7th position of each contig, positions cannot be polymorphic (TGCAGG) which would distort diversity estimates.
# # # So I am going to filter out those positions:
# zcat ParEry.noSEgt2.nogtQ99Cov.noDUST.1.15.depth.gz | sort -Vk1,1 -k2,2 | gawk '$2 !~ /^[234567]$/' > ParEry.noSEgt2.nogtQ99Cov.noDUST.1.15.noTGCAGG.sorted.sites
# 
# 
# # # I want to know how many RAD contigs do these sites belong to:
# cut -f 1 ParEry.noSEgt2.nogtQ99Cov.noDUST.1.15.noTGCAGG.sorted.sites | uniq | wc -l
# # 55100 
# # This looks much better!!!
# 
# #
# # REVAMP OF COVERAGE FILTERING
# #
# # I have written a Python programme that does excess coverage filtering. It filters whole contigs, not single sites.
# # Only SE read coverage at position 2 is counted. So mismapping SE reads and paired end reads are ignored. The pro-
# # gramme determines the individual coverage distributions over contigs that have at least one SE read mapping, i. e.
# # the 0 count category is not included! It then determines the specified percentiles of the individual and across-
# # sample coverage distributions and determines a list of contigs that are below all individual and the across-sample
# # percentile of the coverage distributions. It prints out a new BED file with intervals from contigs that passed the
# # excess coverage filters. Contigs that do not get any SE read mapping at position 2 are also automatically filtered out.
# # Thus, the following command line should be able to replace the whole section on excess coverage filtering above (in
# # addition it also does across sample excess coverage filitering).
# ./excess_coverage_filter.py -b Big_Data_Contigs.noSEgt2.noDUST.bed -p 99 *sorted.bam
# 
# # remove TGCAGG from BED intervals (for description see line 1806)
# perl -lane'if($F[1] == 1){$F[1]=7; print join("	", @F)}else{print}' Big_Data_Contigs.noSEgt2.noDUST.COVfiltered.bed | gawk '$2 < $3' > Big_Data_Contigs.noSEgt2.noDUST.COVfiltered.noTGCAGG.bed
# # note the perl 'join' function needs a literal TAB character
# 
# # create new samtools depth file
# samtools depth  -b Big_Data_Contigs.noSEgt2.noDUST.COVfiltered.noTGCAGG.bed -Q 5 *sorted.bam | bgzip > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.depth.gz
# # Note, that I am skipping '-aa' from the upper command line, so the depth file only contains sites with at least 1x coverage in at least 1 individual.

# # # wrote even coverage filter script in Python
# # cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# # 
# zcat ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.depth.gz | ./minimum_coverage_filter.py -mc 1 -mi 15 | gzip > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sites.gz
# # This takes 2h:35min to finish.
# # The new sites file contains 2,683,395 sites from 85,488,084 sites (3.14%).
# # The filtered sites lie on 52,042 contigs.


# # --------------------
# # HWE filtering
# # --------------------
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# mkdir SNPSTAT
# ln -s ../../BOWTIE2/BAM_dedup/ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sites.gz .
# gzip -dc ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sites.gz | sort -Vk1,1 -k2,2 > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites

# # index sites file for ANGSD
# angsd sites index ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites

# # create regions file for ANGSD
# cut -f 1 ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites | uniq > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.rf

# # split the regions file into small files, each containing only 500 contig names:
# split -l 500 ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.rf SPLIT_RF/

# #
# # estimate total sample F:
# #
# ls SPLIT_RF/* | \
# 	parallel -j 12 "angsd -b bamfile.list -doSnpStat 1 -doMaf 1 -domajorminor 1 -skipTriallelic 1 -gl 1 -snp_pval 1e-2 \
# 	-sites ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SNPSTAT/{/} -rf {}"

# # combine split output files:
# cd SNPSTAT
# zcat aa.hwe.gz | head -n 1 > ParEry.hwe
# for f in *hwe.gz; do zcat $f | tail -n +2 >> ParEry.hwe; done

# #
# # filter for sites with negative F and p-value below 0.05:
# #
# gawk '$7<0 && $9 < 0.05' ParEry.hwe | cut -f 1 | uniq > ParEry.negFis.contigs
# # this finds 439 contigs

# #
# # estimate within population F
# #
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# for POP in PAR ERY; 
# do 
# 	ls SPLIT_RF/* | \
# 		parallel -j 12 "angsd -bam $POP.bamfile.list  -dosnpstat 1 -domaf 1 -domajorminor 1 -skiptriallelic 1 -gl 1 -snp_pval 1e-2 \
# 		-sites ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SNPSTAT/$POP.{/} -rf {}"; 
# done


# # combine split output files:
# cd SnpStat
# for POP in ERY PAR;
# do
# 	zcat $POP.aa.hwe.gz | head -n 1 > $POP.hwe
# 	for f in $POP*hwe.gz; do zcat $f | tail -n +2 >> $POP.hwe; done
# done


# #
# # filter out sites with negative F and p-value below 0.05:
# #
# for POP in ERY PAR;
# do
# 	gawk '$7<0 && $9 < 0.05' $POP.hwe | cut -f 1 | uniq > $POP.negFis.contigs;
# done
# # This finds 173 contigs in PAR and 557 in ERY.

# #
# # filter positive Fis contigs
# #
# # only within pop Fis!
# for POP in ERY PAR;
# do
# 	gawk '$7 > 0 && $9 < 0.05' $POP.hwe | cut -f 1 | uniq > $POP.posFis.contigs;
# done
# # this finds 220 contigs in ERY and 5 in PAR.

# # combine negFis.contigs files:
# cat *.contigs | sort -V | uniq > combined.negPosFis.contigs
# # There are 3486 contigs in the combined list of contigs.


# #
# # create new keep.sites file:
# #
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# perl -ne'chomp; $H{$_}=1; \
# 	END{open(I, "ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites"); while(<I>){@line=split; print if not exists $H{$line[0]}}}' \
# 	SNPSTAT/combined.negPosFis.contigs > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.negPosFisFiltered.sorted.sites
# # I have checked that the new keep.sites file contains 3486 fewer contigs than the one before.
# # The new keep.sites file contains 2,455,851 sites on 48,556 contigs.


# # index sites file for ANGSD
# angsd sites index ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.negPosFisFiltered.sorted.sites
# 
# # create regions file for ANGSD
# cut -f 1 ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.negPosFisFiltered.sorted.sites | uniq > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.negPosFisFiltered.sorted.rf
# 
# 
# # -------------------------------------------------------
# # Sample Allele frequency likelihoods (SAF's) UNFOLDED
# # -------------------------------------------------------
# # # preparation:
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# #mkdir Saf/Ery mkdir Saf/Par
# 
# export keep="ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.negPosFisFiltered.sorted"
# 
# angsd -bam PAR.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out Saf/Par/Par.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1 &
# angsd -bam ERY.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out Saf/Ery/Ery.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1 &
# # note, running with minInd of 9
# # running both commands in parallel (&) takes less than 15 min
# 
# # The SAF file for PAR contains 2,448,252 sites on ... contigs.
# # The SAF file for ERY contains ... sites on ... contigs.
# 
# # --------------------------
# # estimate ML UNFOLDED SFS
# # --------------------------
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# mkdir -p sfs/par sfs/ery
# # 
# realSFS -P 12  Saf/Ery/Ery.unfolded.saf.idx 2>/dev/null > sfs/ery/Ery.unfolded.sfs
# realSFS -P 12  Saf/Par/Par.unfolded.saf.idx 2>/dev/null > sfs/par/Par.unfolded.sfs
# # # this takes less than 15 min to finish
# # # === the spectra are as unusable as with previous filtering ===


# # ===> try with filtering for higher coverage <===
# 
# cd /data3/claudius/Big_Data/BOWTIE2/BAM_dedup
# # 
# zcat ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.depth.gz | ./minimum_coverage_filter.py -mc 3 -mi 10 | gzip > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.depth.gz
# # This takes 2h:35min to finish.
## The Python script minimum_coverage_filter.py should be scraped and replaced by a simple gawk command, which is MUCH faster:
# zcat ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.depth.gz | gawk '{for(i=3; i<=NF; i++) if($i>=3) mi++; if(mi>=10) print; mi=0}' | gzip > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.depth.gz
# # The new sites file contains 414,122 sites from 85,488,084 sites (0.48%). There are only 125,875 sites with 3x coverage in at least 15 individuals.
# # The filtered sites lie on 10,216 contigs.

# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# mkdir SNPstat
# ln -s ../../BOWTIE2/BAM_dedup/ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.sites.gz .
# gzip -dc ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.sites.gz | sort -Vk1,1 -k2,2 > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.sorted.sites

# # index sites file for ANGSD
# angsd sites index ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites


# # create regions file for ANGSD
# cut -f 1 ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.sites | uniq > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.rf

# mkdir SPLIT_rf
# # split the regions file into small files, each containing only 500 contig names:
# split -l 500 ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.1.15.sorted.rf SPLIT_RF/


# #
# # estimate total sample F:
# #
# ls SPLIT_rf/* | \
# 	parallel -j 12 "angsd -b bamfile.list -doSnpStat 1 -doMaf 1 -domajorminor 1 -skipTriallelic 1 -gl 1 -snp_pval 1e-2 \
# 	-sites ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.sorted.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SNPstat/{/} -rf {}"

# 
# # combine split output files:
# cd SNPstat
# zcat aa.hwe.gz | head -n 1 > ParEry.hwe
# for f in *hwe.gz; do zcat $f | tail -n +2 >> ParEry.hwe; done
# 
# #
# # filter for sites with negative F and p-value below 0.05:
# #
# gawk '$7<0 && $9 < 0.05' ParEry.hwe | cut -f 1 | uniq > ParEry.negFis.contigs
# # this finds 259 contigs


# #
# # estimate within population F
# #
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# for POP in PAR ERY; 
# do 
# 	ls SPLIT_rf/* | \
# 		parallel -j 12 "angsd -bam $POP.bamfile.list  -dosnpstat 1 -domaf 1 -domajorminor 1 -skiptriallelic 1 -gl 1 -snp_pval 1e-2 \
# 		-sites ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.sorted.sites -only_proper_pairs 0 -minMapQ 5 -minQ 20 -out SNPstat/$POP.{/} -rf {}"; 
# done


# #
# # filter out sites with negative F and p-value below 0.05:
# #
# for POP in ERY PAR;
# do
# 	gawk '$7<0 && $9 < 0.05' $POP.hwe | cut -f 1 | uniq > $POP.negFis.contigs;
# done
# # This finds 113 contigs in PAR and 322 in ERY.


# #
# # filter positive Fis contigs
# #
# # only within pop Fis! We don't want to filter out divergent SNP's.
# for POP in ERY PAR;
# do
# 	gawk '$7 > 0 && $9 < 0.05' $POP.hwe | cut -f 1 | uniq > $POP.posFis.contigs;
# done
# # this finds 346 contigs in ERY and 326 in PAR.


# # combine Fis.contigs files:
# cat *Fis.contigs | sort -V | uniq > combined.negPosFis.contigs
# # There are 1097 contigs in the combined list of contigs.


# # #
# # # create new keep.sites file:
# # #
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# perl -ne'chomp; $H{$_}=1; \
# 	END{open(I, "ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.sorted.sites"); while(<I>){@line=split; print if not exists $H{$line[0]}}}' \
# 	SNPstat/combined.negPosFis.contigs > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.negPosFisFiltered.sorted.sites
# # I have checked that the new keep.sites file contains 1097 fewer contigs than the one before.
# # The new keep.sites file contains 368,764 sites on 9,119 contigs.

# # index sites file for ANGSD
# angsd sites index ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.negPosFisFiltered.sorted.sites


# # create regions file for ANGSD
# cut -f 1 ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.negPosFisFiltered.sorted.sites | uniq > ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.negPosFisFiltered.sorted.rf


# # -------------------------------------------------------
# # Sample Allele frequency likelihoods (SAF's) UNFOLDED
# # -------------------------------------------------------
# # # preparation:
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# # #mkdir -p saf/ery saf/par
# # 
# export keep="ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.negPosFisFiltered.sorted"
# 
# angsd -bam PAR.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out saf/par/par.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1 2>/dev/null &
# angsd -bam ERY.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out saf/ery/ery.unfolded -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1 2>/dev/null &
# # note, running with minInd of 9

# # --------------------------
# # estimate ML UNFOLDED SFS
# # --------------------------
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# mkdir -p afs/par afs/ery
# realSFS -P 10 -maxIter 50000 -tole 1e-6 -m 0 saf/ery/ery.unfolded.saf.idx 2>/dev/null > afs/ery/ery.unfolded.sfs &
# realSFS -P 10 -maxIter 50000 -tole 1e-6 -m 0 saf/par/par.unfolded.saf.idx 2>/dev/null > afs/par/par.unfolded.sfs &


# Let's enfore at least 15 individuals with read data:

# # -------------------------------------------------------
# # Sample Allele frequency likelihoods (SAF's) UNFOLDED
# # -------------------------------------------------------
# # # preparation:
# cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
# export keep="ParEry.noSEgt2.noDUST.COVfiltered.noTGCAGG.3.10.negPosFisFiltered.sorted"
# 
# angsd -bam PAR.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out saf/par/par.unfolded.15 -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 15 -GL 1 -doSaf 1 -nThreads 1 2>/dev/null &
# angsd -bam ERY.bamfile.list -ref Big_Data_ref.fa -anc Big_Data_ref.fa -out saf/ery/ery.unfolded.15 -fold 0 \
# 	-sites $keep.sites -rf $keep.rf -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 15 -GL 1 -doSaf 1 -nThreads 1 2>/dev/null &
# # note, running with minInd of 9

# # --------------------------
# # estimate ML UNFOLDED SFS
# # --------------------------
cd /data3/claudius/Big_Data/ANGSD/DEDUPLICATED
realSFS -P 10 -maxIter 50000 -tole 1e-6 -m 0 saf/ery/ery.unfolded.15.saf.idx 2>/dev/null > afs/ery/ery.unfolded.15.sfs &
realSFS -P 10 -maxIter 50000 -tole 1e-6 -m 0 saf/par/par.unfolded.15.saf.idx 2>/dev/null > afs/par/par.unfolded.15.sfs &
# # ===> END with filtering for higher coverage <===




