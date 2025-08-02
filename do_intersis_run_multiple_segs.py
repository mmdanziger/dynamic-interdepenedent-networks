try:
    import intersis
except ImportError:
    from oscillators import intersis
from sys import argv
import numpy as np
import json
from os import environ


def num(x):
    try:
        return int(x)
    except:
        return float(x)

N = int(argv[1])
k = num(argv[2])
f = num(argv[3])
interaction_idx = int(argv[4])
beta = np.array([num(argv[5]),num(argv[6])])
thetamin,thetamax = [num(argv[7]),num(argv[8])]
steps = int(argv[9])
time_per_seg = num(argv[10])
use_local_state = environ.get("LOCAL_STATE")
sis = intersis.InterSIS(N,k,beta,"er")

if use_local_state:
    sis.local_mf = False
sis.interaction_type = ["interdependent","competitive","hybrid"][interaction_idx]
sis.hybrid = True if sis.interaction_type == "hybrid" else False


res = {}
for theta1 in np.linspace(thetamin,thetamax,steps):
    for theta2 in np.linspace(thetamin, thetamax, steps):
        sis.randinit_upto((theta1*2,theta2*2))
        segs = list(map(lambda x: [list(i) for i in x],sis.integrated_segments(6,time_per_seg)))
        this_key = "(%.4f,%.4f)"%(theta1,theta2)
        res[this_key]=segs
        print(this_key)
if use_local_state:
    fname_base = "beta_samplesLS_"
else:
    fname_base = "beta_samples_"
json.dump(res,open(fname_base + "_".join(argv[1:])+".json","w"))
