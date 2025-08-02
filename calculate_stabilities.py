from __future__ import division
import single_net_sync as sns
import numpy as np
import json
from sys import argv
import progressbar

def get_pbar(frac_str_rval):
    widgets = ["All f  ", progressbar.Percentage(),
               ' ', progressbar.Bar(),
               ' ', progressbar.ETA(),
               ' ', progressbar.AdaptiveETA()]
    oldwidgets = ['Test: ', progressbar.Percentage(), ' ', progressbar.Bar(marker=progressbar.RotatingMarker()),
                  ' ', progressbar.ETA(), ' ', progressbar.FileTransferSpeed()]
    return progressbar.ProgressBar(maxval=len(frac_str_rval), widgets=widgets).start()


k = 50
lam = 0.14
frac = 0.1
eps = 1e-12


def get_stability_vec(frac_str_r_val_list):
    global solutions,k
    frac_str,r_val_list = frac_str_r_val_list
    stabilities=[]
    for (lam_val, r_vals) in zip(lambda_vec, r_val_list):
        stable_vec = list(map(lambda r: sns.get_energy_level(intfunc, lam_val, float(frac_str), k, r), r_vals))
        stabilities.append(stable_vec)
    print("f = %s complete"%frac_str)
    return [frac_str,stabilities]

def my_callback(x):
    global solutions
    solutions["energy_level"] = dict( (frac_str,stabilities) for frac_str,stabilities in x)

intfunc = sns.get_interpolated_integral(k, mode="r",fname="interpolation_values_combined_1d_k50_HD.json")
results_fname = argv[1] if len(argv) > 1 else "zhang_results_k50_HD_superzoom.json"
solutions = json.load(open(results_fname))
lambda_vec = solutions["lambda_vec"]
solutions["energy_level"] = {}

frac_str_r_val_list = sorted(solutions["r"].items(),key=lambda x: float(x[0]))
from multiprocessing import Pool
mypool = Pool(processes=8)
#bar = get_pbar(frac_str_r_val_list)
r = mypool.map_async(get_stability_vec,frac_str_r_val_list,callback=my_callback)
bar_state=0
r.wait()

# for frac_str, r_val_list in sorted(solutions["r"].items(),key=lambda x: float(x[0])):
#     solutions["energy_level"][frac_str] = []
#     bar = get_pbar(frac_str,lambda_vec)
#     lam_r_vals = list(zip(lambda_vec, r_val_list))
#     for idx, (lam_val, r_vals) in enumerate(zip(lambda_vec, r_val_list)):
#         stable_vec = list(map(lambda r: sns.get_energy_level(intfunc, lam_val, float(frac_str), 12, r), r_vals))
#         solutions["energy_level"][frac_str].append(stable_vec)
#         bar.update(idx)
#     bar.finish()

json.dump(solutions, open("zhang_results_energy_k50_HD_superzoom.json", "w"))
