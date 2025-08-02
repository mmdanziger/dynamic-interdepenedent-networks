from __future__ import division,print_function
import two_nets_sync as tns
from sys import argv
from collections import defaultdict
import json
from glob import glob
import os

"""
Script to generate solutions for a given value of f.
It can run as either the master or the slave.
The master has no args, and launches slaves with args indicating which values to cover.
This way multiprocessing can be used cleanly and efficiently.
After all of the slaves write their output to /tmp, the master compiles all the outputs into one json.
"""

k = 12
frac = 1
lambda_step = 0.001
lambda_i = 0
lambda_f = 0.5
interaction_type="competitive"
if len(argv) > 1:
    """
    Slave
    """
    lambda1_i = float(argv[1])
    lambda1_f = float(argv[2])
    res = defaultdict(dict)
    integral_fname = "/home/micha/phd/oscillators/interpolation_values_combined_1d_k%i_HD.json" % k
    try:
        os.stat(integral_fname)
    except:
        integral_fname = "/storage" + integral_fname

    intfunc = tns.sns.get_interpolated_integral(kbar=k, mode="r",fname=integral_fname)
    for lambda1 in tns.np.arange(lambda1_i, lambda1_f, step=lambda_step):
        for lambda2 in tns.np.arange(lambda_i, lambda_f, step=lambda_step):
            r1sols, r2sols = tns.solve_R1R2(intfunc, [lambda1, lambda2], frac,interaction_type=interaction_type)
            res["%.8f" % lambda1]["%.8f" % lambda2] = {"r1": list(r1sols), "r2": list(r2sols)}
    json.dump(res,
              open("/tmp/two_net_result_%s_k%i_f%.4f_lambdai%.8f_lambdaf%.8f.json" % (interaction_type, k, frac, lambda1_i, lambda1_f), "w"))
else:
    """
    Master
    """
    from subprocess import call

    def run(li, lf):
        call(["python", argv[0], str(li), str(lf)])
        print("complete.")


    import multiprocessing

    processes = multiprocessing.cpu_count()
    jobs = []
    for i in range(processes):
        li = lambda_i + i * (lambda_f - lambda_i) / processes
        lf = lambda_i + (i + 1) * (lambda_f - lambda_i) / processes
        p = multiprocessing.Process(target=run, args=(li, lf))
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()
    """
    Collect and save data from slaves' outputs.
    """
    d = {}
    for fname in glob("/tmp/two_net_result_%s_k%i_f%.4f_*.json" % (interaction_type,k, frac)):
        subd = json.load(open(fname))
        for key,v in subd.items():
            d[key] = v
    json.dump(d,open("two_net_result_%s_k%i_f%.4f.json"%(interaction_type,k,frac),"w"))
