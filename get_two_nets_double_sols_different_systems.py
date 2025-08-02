from __future__ import division,print_function
try:
    import two_nets_sync as tns
except ImportError:
    from oscillators import two_nets_sync as tns
from sys import argv
from collections import defaultdict
import json
from glob import glob
import os
from os.path import basename

"""
Script to generate solutions for a Kuramoto-SIS mixed system and a given value of f.
It can run as either the master or the slave.
The master has no args, and launches slaves with args indicating which values to cover.
This way multiprocessing can be used cleanly and efficiently.
After all of the slaves write their output to /tmp, the master compiles all the outputs into one json.
"""

def num(string_input):
    try:
        return int(string_input)
    except:
        return float(string_input)
force_fname = False
mode = argv[1]
try:
    k1 = num(argv[2])#or gamma if sf
    k2 = num(argv[3])
except:
    integral1_fname = argv[2]
    integral2_fname = argv[3]
    force_fname = True
    k1 = [num(i[1:]) for i in basename(integral1_fname).split("_") if i[0]=="k"][0]
    k2 = [num(i[1:]) for i in basename(integral2_fname).split("_") if i[0] == "k"][0]

frac = float(argv[4])
interaction_idx = int(argv[5])
lambda_i = float(argv[6])
lambda_f = float(argv[7])
lambda_step = float(argv[8])
degree_dist = "er"
q = float(os.environ["OSCILLATORS_Q"]) if "OSCILLATORS_Q" in os.environ else 0

interaction_type=["interdependent","competitive", "hybrid","mixed","halfhalf"][interaction_idx]

if mode == "slave":
    """
    Slave
    """
    lambda1_i = lambda_i
    lambda1_f = lambda_f
    lambda2_i = float(argv[9])
    lambda2_f = float(argv[10])
    res = defaultdict(dict)
    if force_fname:
        integral1_fname = integral1_fname
        integral2_fname = integral2_fname
    elif degree_dist == "er":
        integral1_fname = "/home/micha/phd/oscillators/interpolation_values_combined_1d_k%i_HD.json" % k1
        integral2_fname = "/home/micha/phd/oscillators/sis_interpolation_values_combined_1d_k%i_HD.json" %k2
    elif degree_dist == "sf":
        k1_string = "%i"%k1 if type(k1) is int else "%.1f"%k1
        integral1_fname = "/home/micha/phd/oscillators/interpolation_values_combined_1d_gamma%s_HD.json" % k1_string
    try:
        os.stat(integral1_fname)
        os.stat(integral2_fname)

    except:
        integral1_fname = "/storage" + integral1_fname
        integral2_fname = "/storage" + integral2_fname

    intfunc1 = tns.sns.get_interpolated_integral(kbar=k1, mode="r",fname=integral1_fname,gamma=k1)
    intfunc2 = tns.sns.get_interpolated_integral(kbar=k2, mode="r", fname=integral2_fname, gamma=k2)

    for lambda1 in tns.np.arange(lambda1_i, lambda1_f, step=lambda_step):
        for lambda2 in tns.np.arange(lambda2_i, lambda2_f, step=lambda_step):
            sols = tns.solve_R1R2_different_systems_double_solutions(intfunc1, intfunc2, [lambda1, lambda2], frac,interaction_type=interaction_type,q=q)
            res["%.8f" % lambda1]["%.8f" % lambda2] = sols
    json.dump(res,
              open("/tmp/two_net_different_systems_result_%s_k1%i_k2%i_f%.4f_lambdai%.8f_lambdaf%.8f.json" % (interaction_type, k1, k2, frac, lambda1_i, lambda1_f), "w"))
else:
    """
    Master
    """
    from subprocess import call

    cores = int(argv[9]) if len(argv) > 9 else None


    def run(li, lf):
        global force_fname
        k1_arg = integral1_fname if force_fname else str(k1)
        k2_arg = integral2_fname if force_fname else str(k2)

        call(["python", argv[0], "slave", k1_arg, k2_arg, "%.4f"%frac, str(interaction_idx), str(li), str(lf),str(lambda_step),str(lambda_i),str(lambda_f)])
        print("complete.")


    import multiprocessing

    processes = multiprocessing.cpu_count() if not cores else cores
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
    for fname in glob("/tmp/two_net_different_systems_result_%s_k1%i_k2%i_f%.4f_*.json" % (interaction_type,k1,k2, frac)):
        subd = json.load(open(fname))
        for key,v in subd.items():
            d[key] = v
        call(["rm" , fname])

    if degree_dist == "er":
        final_fname = "two_net_different_systems_result_%s_k1%i_k2%i_f%.4f_q%.4f.json" % (interaction_type, k1, k2, frac, q)
    else:
        final_fname = "two_net_different_systems_result_%s_gamma1%i_gamma2%i_f%.4f_q%.4f.json" % (interaction_type, k1,k2, frac,q)

    json.dump(d,open(final_fname,"w"))
    try:
        call(["scp", final_fname, "micha@thedanziger.com:/home/micha"])
    except:
        pass
