#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: print_Q99_exCov_contigs.pl
#
#        USAGE: ./print_Q99_exCov_contigs.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 15/11/16 16:30:18
#     REVISION: ---
#===============================================================================

use strict;
use warnings;


use strict;
use warnings;

open(my $Q99, "<Q99.cov") or die $!;

# read in list of 99% coverage percentiles per ind BAM file
my @line = ();
my %fn_Q99 = ();
while(<$Q99>){
	next if /#/;
	@line = split;
	$fn_Q99{$line[0]} = $line[1]; # the hash contains file names as keys and 99% coverage percentiles as value
}
close($Q99);


my %exCov_contigs = ();
for my $fn (keys %fn_Q99){
	# opens a stream of SE reads, cuts out contig name 
	# and collapses BAM records with identical mapping position while keeping the count
	open(my $BAM, "samtools view -f64 $fn | cut -f3 | uniq -c |") or die $!; 	
	while(<$BAM>){
		@line = split;
		if($line[0] > $fn_Q99{$fn}){ # if the uniq count (i. e. coverage) is greater than Q99 for that individual
			$exCov_contigs{$line[1]} = 1; # store contig name 
		}
	}
	close($BAM);
}

open(my $OUT, ">", "Big_Data_Contigs.gtQ99Cov") or die $!;

for my $contig (keys %exCov_contigs){
	print $OUT $contig, "\n";
}
