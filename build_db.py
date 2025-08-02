from __future__ import division,print_function

import lambdac_results_loader as lrl
from sys import argv
from glob import glob

lamc_glob_string = argv[1]
energy_flist_fname = argv[2]

lrl.build_db(lambda_c_glob=lamc_glob_string,trhist_glob="/storage/home/micha/cudamoto_results/trhist/processed/*",energy_fname=None)

for row in open(energy_flist_fname):
    k,fname = row.strip().split()
    lrl.add_sols_to_db(k=int(k),energy_fname=fname)
