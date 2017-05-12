#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: subtract.pl
#
#        USAGE: ./subtract.pl  
#
#  DESCRIPTION: filter BED file
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Claudius Kerth (CEK), c.kerth[at]sheffield.ac.uk
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 15/11/16 18:46:38
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

# this bed file contains intervals spanning the whole Big Data reference
open(my $ALL, "Big_Data_ref.fa.bed") or die $!;

my %ALL = ();
my ($contig, $start, $end) = ("", 0, 0);

# store all intervals in a hash
while(<$ALL>){
	($contig, $start, $end) = split;	
#	print $contig, ":", $start, "-", $end, "\n";
	$ALL{$contig} = "\t$start\t$end\n"; # key: contig name; value: <TAB>interval_start<TAB>interval_end<LF>
}
close($ALL);
#for $contig (keys %ALL){
#	print $contig, $ALL{$contig};
#}

# # this is a list of unique contig names that occured in the BAM files
# open(my $WM, "Big_Data_Contigs_with_mappings.list") or die $!;
# 
# my %WM = ();
# 
# # store contig names in a new hash and provide interval coordinates as values
# while(<$WM>){
# 	chomp;
# 	$WM{$_} = $ALL{$_}; 
# }
# #foreach $contig (keys %WM){
# #	print $contig, $WM{$contig};
# #}
# undef %ALL; # clear memory
# 
# #print scalar keys %WM, "\n";

# this file contains a list of contig names with SE read mapping beyond position 2
# they should be removed from the %WM hash
open(my $mismap, "Big_Data_Contigs.SEposgt2") or die $!;

while(<$mismap>){
	chomp;
	delete $ALL{$_};
}
close($mismap);

#print scalar keys %WM, "\n";

# # this contains a list of contig names with coverage above 99th percentile
# # they should be removed from %WM
# open(my $exCov, "Big_Data_Contigs.gtQ99Cov") or die $!;
# 
# while(<$exCov>){
# 	chomp;
# 	delete $WM{$_};
# }

#print scalar keys %WM, "\n";

# print the remaining BED intervals to STDOUT
foreach $contig (keys %ALL){
	print $contig, $ALL{$contig};
}












