try:
    import intersis
except ImportError:
    from oscillators import intersis
from sys import argv
import numpy as np
from os import environ
from datetime import datetime

def num(x):
    try:
        return int(x)
    except:
        return float(x)

N = int(argv[1])
k = num(argv[2])
f = num(argv[3])
interaction_idx = int(argv[4])
p0 = np.array([num(argv[5]),num(argv[6])])
p1 = np.array([num(argv[7]),num(argv[8])])
steps = int(argv[9])
tfinal = num(argv[10])
timestamp = datetime.now().isoformat()
use_local_state = environ.get("LOCAL_STATE")
use_global_op = environ.get("GLOBAL_OP")
sis = intersis.InterSIS(N,k,[1,1],"er")
if use_local_state:
    sis.local_mf = False

sis.hybrid = False
sis.interaction_type = ["interdependent","competitive","hybrid"][interaction_idx]
sis.hybrid = True if sis.interaction_type == "hybrid" else False
sis.randinit_upto([1,1])
sis.set_f(f)
if use_global_op:
    bv, hv, ghv = sis.betapath_from_points_global_op_compare(p0, p1, steps, tfinal)
else:
    bv, hv = sis.betapath_from_points(p0, p1, steps, tfinal)

fname_base = "betapath"
if use_local_state:
    fname_base += "LS"

np.save(fname_base + "_" + "_".join(argv[1:] + [timestamp])+".npy", np.array([bv,hv]))

if use_global_op:
    fname_base += "GLOBALOP"
    np.save(fname_base + "_" + "_".join(argv[1:] + [timestamp]) + ".npy", np.array([bv, ghv]))

