#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: slim.BAMs.pl
#
#        USAGE: ./slim.BAMs.pl  
#
#  DESCRIPTION: creates new BAM files only containing SAM records for filtered contigs
#               this was an attempt to speed up ANGSD's region lookup
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 16/11/16 11:33:03
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

# open region file made for ANGSD that contains all contigs with sites to analyse:
open(my $CONTIGS, "../ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.noNegFis.sorted.rf") or die $!;

my %KeepContigs = ();

# store contig names in hash table:
while(<$CONTIGS>){
	chomp;
	s/://;
	$KeepContigs{$_} = 1;
}
close($CONTIGS);


# now go through every BAM file and print out SQ header lines and BAM records
# that belong to the above contigs

my $contig = "";
my @line = ();

foreach my $bam (glob("*sorted.slim.bam")){
	# open stream of uncompressed SAM
	open( my $BAM, "samtools view -h $bam | ") or die $!;
	$bam =~ s/\.slim\.bam$//;
	# write output to samtools view for conversion to BAM
	open( my $SLIM, "| samtools view -b - > $bam.SLIM.bam" ) or die $!;

	while(<$BAM>){
		# if SAM header line
		if(/^\@/){
			# if ref seq line
			if(/^\@SQ/){
				# get contig name
				($contig) = $_ =~ m/SN:(Contig_\d+)\t/;
				# print header line if contig name is in the previously created hash table
				print $SLIM $_ if exists $KeepContigs{$contig};
			}
			# if any other header line
			else{ print $SLIM $_; }
		}
		# if SAM record
		else{
			# split line into 4 parts
			@line = split("\t", $_, 4);
			# print SAM record if the contig is in the hash table
			print $SLIM $_ if exists $KeepContigs{$line[2]};
		}
	}
	close($BAM);
	close($SLIM);
}
