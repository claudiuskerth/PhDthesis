def flatten(array):
    """
    Returns a list of flattened elements of every inner lists (or tuples)
    ****RECURSIVE****
    """
    import numpy
    res = []
    for el in array:
        if isinstance(el, (list, tuple, numpy.ndarray)):
            res.extend(flatten(el))
            continue
        res.append(el)
    return list( res )


# ====================================================


def get_flag_count(out, NM=True):
    """
    out: list of tuples, each containing p_init and popt + additional info, including warnflags
    as produced by run_dadi.py
    """

    from collections import defaultdict

    if NM: # if ar from Nelder-Mead
        i = 4 # the warnflag is reported at index position 4 in the output array
    else: # ar from BFGS optimisation
        i = 6

    warnflag = defaultdict(int)

    for res in out:
        if res[1][i] == 1: # notice the change in indexing
            warnflag[1] += 1
        elif res[1][i] == 2:
            warnflag[2] += 1
        elif res[1][i] == 0:
            warnflag[0] += 1
        else:
            warnflag[999] +=1

    if NM:
        print "success", warnflag[0]
        print "Maximum number of function evaluations made.", warnflag[1]
        print "Maximum number of iterations reached.", warnflag[2]
        print "unknown flag", warnflag[999]
    else:
        print "success", warnflag[0]
        print "Maximum number of iterations exceeded.", warnflag[1]
        print "Gradient and/or function calls not changing.", warnflag[2]
        print "unknown flag", warnflag[999]

# ====================================================
# 
# # import global namespace from calling script to access its global names:
# from __main__ import *
# # this is because this module has its own global namespace
# 
# def run_dadi(p_init): # for the function to be called with map, it needs to have one input variable
#     """
#     p_init: initial parameter values to run optimisation from
#     """
# 
#     if perturb == True:
#         p_init = dadi.Misc.perturb_params(p_init, fold=fold, upper_bound=upper_bound, lower_bound=lower_bound)
#         # note upper_bound and lower_bound variables are expected to be in the namespace of each engine
# 
#     # run optimisation of paramters
#     popt = dadi_opt_func(p0=p_init, data=sfs, model_func=func_ex, pts=pts_l, \
#     lower_bound=lower_bound, upper_bound=upper_bound, \
#     verbose=verbose, maxiter=maxiter, full_output=full_output)
# 
#     # pickle to file
#     import dill
#     name = outname[:] # make copy of file name stub!
#     for p in p_init:
#         name += "_%.4f" % (p)
#     with open(name + ".dill", "w") as fh:
#     dill.dump((p_init, popt), fh)
# 
#     return p_init, popt
