from __future__ import division,print_function
import single_net_sync as sns
import numpy as np
import json

import progressbar

def get_pbar(frac_str,lambda_vec):
    widgets = ["f = " + frac_str, progressbar.Percentage(),
               ' ', progressbar.Bar(),
               ' ', progressbar.ETA(),
               ' ', progressbar.AdaptiveETA()]
    oldwidgets = ['Test: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()),
                  ' ', progressbar.ETA(), ' ', progressbar.FileTransferSpeed()]
    return progressbar.ProgressBar(maxval=len(lambda_vec), widgets=widgets).start()


def get_sols(f):
    global k,lambda_vec,intfunc
    sols = []
    for idx,lam in enumerate(lambda_vec):
        sols.append(sns.solve_R(intfunc, lam, frac, k))
    print("Completed f = %.4f"%f)
    return sols


def callback(sols):
    global fvec,results
    for frac,sol in zip(fvec,sols):
        results["%.4f"%frac] = sol

k = 50
lambda_vec = np.linspace(0.0200111, 0.03488555, 2000)

from multiprocessing import Pool


lam = 0.14
frac = 0.1
eps = 1e-12
print("Interpolation...", end='')
intfunc = sns.get_interpolated_integral(k,mode="r",fname="interpolation_values_combined_1d_k50_HD.json")
print("complete.")
results = {}
fvec = [0.01 * i for i in range(101)]
mypool = Pool(processes=8)

r = mypool.map_async(get_sols, fvec,callback = callback)

r.wait()

# for frac in fvec:
#     bar = get_pbar(str(frac), lambda_vec)
#     results["%.4f"%frac] = []
#     for idx,lam in enumerate(lambda_vec):
#         results["%.4f"%frac].append(sns.solve_R(intfunc, lam, frac, 12))
#         bar.update(idx)
#     bar.finish()
print("complete.")
json.dump({"lambda_vec": [i for i in lambda_vec], "k" : k , "r": results}, open("zhang_results_k50_HD_superzoom.json", "w"))
