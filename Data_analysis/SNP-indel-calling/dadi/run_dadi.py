#!/usr/bin/env python

# ----------------------------------------
# is this using the right python version?
# ----------------------------------------
# import sys
# print sys.version_info

# ----------------------------------------
# import modules
# ----------------------------------------
import sys

sys.path.insert(0, '/home/claudius/Downloads/dadi')

import dadi
import numpy


# ----------------------------------------
# read in command line arguments
# ----------------------------------------
import argparse

parser = argparse.ArgumentParser(description="runs dadi optimisation of parameters for given model")
parser.add_argument("-p", "--path_to_spectrum_file", help="file path to site frequency data")
parser.add_argument("-m", "--dadi_model", help="model function to use")
parser.add_argument("-u", "--upper", help="upper parameter bound")
parser.add_argument("-l", "--lower", help="lower parameter bound")
parser.add_argument("-i", "--p_init", help="initial parameter values")
parser.add_argument("-d", "--dadi_opt_func", help="dadi optimisation function to use", default="dadi.Inference.optimize_log")
parser.add_argument("-s", "--stub", help="file name stub for output files", default="out_run_dadi")
parser.add_argument("--maxiter", help="maximum number of iterations allowed", default=10, type=int)
args = parser.parse_args()

# ----------------------------------------
# load spectrum from file
# ----------------------------------------
# check that file exists at the path
import os
assert os.path.exists(args.path_to_spectrum_file), "spectrum file cannot be found at the given location"
sfs = dadi.Spectrum.from_file(args.path_to_spectrum_file)

# ----------------------------------------
# set up 
# ----------------------------------------

ns = sfs.sample_sizes

# grid sizes
pts_l = [ns+10, ns+20, ns+30]

func = eval(args.dadi_model)

func_ex = dadi.Numerics.make_extrap_log_func(func)

upper_bound = eval(args.upper)
lower_bound = eval(args.lower)

# print upper_bound
# print lower_bound


# ----------------------------------------
# optimisation
# ----------------------------------------
dadi_opt_func = eval(args.dadi_opt_func)
p_init = eval(args.p_init)
popt = dadi_opt_func(p0=p_init, data=sfs, model_func=func_ex, pts=pts_l, \
        lower_bound=lower_bound, upper_bound=upper_bound, \
        verbose=0, maxiter=args.maxiter, full_output=True)

# # ----------------------------------------
# # pickle optimal parameters to file
# # ----------------------------------------
import pickle

# print p_init
outname = args.stub
for p in p_init:
    if p < 1:
        outname += "_%f" % (p)
    else:
        outname += "_%i" % (p)
outname += ".pickle"
# print outname

fh = open("OUT_run_dadi/" + outname, "w")
pickle.dump((p_init, popt), fh)
fh.close()

