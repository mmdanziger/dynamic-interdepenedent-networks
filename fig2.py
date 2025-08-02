from __future__ import division, print_function
import lambdac_results_loader as lrl
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import json
from os.path import basename

def num(x):
    try:
        return int(x)
    except ValueError:
        return float(x)
def metadata(fname):
    return dict( (i[0],num(i[1:])) for i in basename(fname).strip().split("_") if not i[0].isnumeric() and i[1].isnumeric())
l = lrl.Loader()
simulations_file = argv[1]

flist = sorted(list(map(lambda x: metadata(x)["f"],open(simulations_file).readlines())))
klist = sorted(list(map(lambda x: metadata(x)["k"],open(simulations_file).readlines())))
assert(all([i == klist[0] for i in klist]))
k = klist[0]

f = plt.figure()
l.plot_r_of_lambda(k=k, flist=sorted(list(set(flist))))
thiscolorlist = dict( ("%.2f"%frac,col) for frac,col in zip(sorted(set(flist)),reversed(l.color_list)))
minlam,maxlam = [1,0]
for fname in open(simulations_file):
    c = thiscolorlist["%.2f"%metadata(fname)["f"]]
    print((metadata(fname)["f"],c))
    lamvec,rvec = zip(*(json.load(open(fname.strip()))))
    maxlam = max( maxlam, max(lamvec))
    minlam = min( minlam, min(lamvec))
    plt.plot(lamvec,rvec, 's-', ms=3.5, lw=0.8,color = c, mew=0,alpha=0.5)

plt.xlim([minlam,maxlam])
plt.tight_layout()