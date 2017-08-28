#!/usr/bin/env python

from glob import glob
import re
from collections import Counter
import subprocess32 as sp
import string
from itertools import product
from sys import stderr
from time import time

def split_regions_file(boot_contigs_dict, fnames, size):
    """
    takes Counter dictionary of bootstrapped contigs
    and an iterator over filenames to choose

    writes out split regions files with repetitions of contigs
    NOT spread over different split regions files
    """
    c = 0 # initialise contig count
    # get next file name from iterator
    fn = fnames.next()
    # open new file for writing and get filehandle
    out = open("split_rf/" + fn[0] + fn[1], "w")
    # iterate over Counter dict of bootstrapped contigs, key=contig name, value=count (rep)
    for contig,rep in sorted(boot_contigs_dict.items(), key=lambda x: int(x[0].replace("Contig_", ""))):
        c+=rep
        if c > size: # write up to 'size' contigs to each split rf file
            out.close() # close current rf file
            fn = fnames.next() # get next file name from iterator
            out = open("split_rf/" + fn[0] + fn[1], "w") # open new rf file for writing
            c = rep
        for _ in range(rep): # write contig name to rf file as often as it occurs in the bootstrap resample
            out.write(contig + "\n")


index = '' # index of bootstrap replicate

for rf in sorted(glob("including_non-overlapping/BOOT_RF/000*")):
    start = time()
    index = re.findall(r'\d+', rf)[-1]
    # reset array for bootstrapped contigs
    boot_contigs = [] 
    with open(rf, "r") as boot_rf:
        for contig in boot_rf:
            boot_contigs.append(contig.rstrip())
    # create dictionary of counts of contigs
    boot_contigs_dict = Counter(boot_contigs)
    # clear directory
    sp.call("rm -f split_rf/*", shell=True)
    # get filename iterator
    fnames = product(string.lowercase, repeat=2)
    # split bootstrapped regions file, 400 contigs per file
    split_regions_file(boot_contigs_dict, fnames, 400)
    # remove previous split SAF files for PAR
    cmd = "rm -f including_non-overlapping/SAF/bootstrap/PAR/[a-z]*"
    sp.call(cmd, shell=True)
    # remove previous split SAF files for ERY
    cmd = cmd.replace("PAR", "ERY")
    sp.call(cmd, shell=True)
    # run SAF calculation in parallel for PAR
    cmd = 'ls split_rf/* | parallel -j 24 "angsd -bam PAR.slim.bamfile.list -ref Big_Data_ref.fa \
            -anc Big_Data_ref.fa -out including_non-overlapping/SAF/bootstrap/PAR/{/}.unfolded -fold 0 \
            -sites all.sites -rf {} -only_proper_pairs 0 -baq 1 -minMapQ 5 -minInd 9 -GL 1 -doSaf 1 -nThreads 1 2>/dev/null"'
    sp.call(cmd, shell=True)
    # run SAF calculation in parallel for ERY
    cmd = cmd.replace("PAR", "ERY")
    sp.call(cmd, shell=True)
    # concatenate split SAF files for PAR
    cmd = "realSFS cat -outnames including_non-overlapping/SAF/bootstrap/PAR/{}.unfolded including_non-overlapping/SAF/bootstrap/PAR/[a-z]*saf.idx 2>/dev/null".format(index)
    sp.call(cmd, shell=True)
    # concatenate split SAF files for ERY
    cmd = cmd.replace("PAR", "ERY")
    sp.call(cmd, shell=True)
    end = time()
    run_time = end - start
    print >> stderr, "Finished SAF calculation for bootstrap {0}. It took {1} sec to complete.".format(index, int(run_time))


# remove split SAF files for PAR
cmd = "rm -f including_non-overlapping/SAF/bootstrap/PAR/[a-z]*"
sp.call(cmd, shell=True)
# remove split SAF files for ERY
cmd = cmd.replace("PAR", "ERY")
sp.call(cmd, shell=True)
