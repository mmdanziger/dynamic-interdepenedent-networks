from __future__ import division,print_function
from oscillators.fig2_8 import get_lifetime
from sys import argv
from glob import glob
from os.path import join
from pyNoN import logbin
import json

dirname = argv[1]
filetype = argv[2]
sis=False
if filetype == "hdf5":
    sis = False
elif filetype == "json":
    sis = True

betac_vector = []
lifetime_list = []

for fname in glob(join(dirname,"*."+filetype)):
    try:
        betac,betalist,lifetime = get_lifetime(fname)
        lifetime_list.append(lifetime)
        betac_vector.append(betac)
    except:
        print("Unable to read data in %s , skipping"%fname)


points = []
betac = sum(map(float,betac_vector)) / len(betac_vector)

for betac_,lifetime in zip(betac_vector,lifetime_list):
    for beta,this_lifetime in lifetime:
        points.append((float(betac_) - float(beta), this_lifetime))

lb = logbin.LogBin(points=points)
ofname = "plateau_average_""_".join(argv[1:])+".json"
json.dump({"xavg":lb.xavg,"yavg":lb.yavg,"xerr":lb.xerr,"yerr":lb.yerr},open(ofname,"w"))

print("beta_c = %.7f"%betac)

