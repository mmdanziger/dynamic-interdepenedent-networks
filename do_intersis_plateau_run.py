try:
    import intersis
except ImportError:
    from oscillators import intersis
from sys import argv
import numpy as np
import json
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
beta0 = num(argv[5])
beta1 = num(argv[6])
deltabeta = num(argv[7])
tfinal = num(argv[8])

sis = intersis.InterSIS(N,k,[1,1],"er")
sis.hybrid = False
sis.interaction_type = ["interdependent","competitive","hybrid"][interaction_idx]
if sis.interaction_type == "hybrid":
    sis.hybrid = True
sis.randinit_upto([1,1])
sis.set_f(f)
beta_min = beta0
beta_max = beta1
outlist=[]
betalist=[]
threshold = 0.001

betac = 4 / k


if False:
    for beta in betalist:
        sis.beta = np.array([beta,beta])
        th1th2 = sis.integrate_by_chunks(tfinal,stop_condition=True)
        outlist.append(list(np.vsplit(th1th2,2)[0].flatten()))
        final = outlist[-1][-1]
        print("%.7f -> %.8f"%(beta,final))
if True:
    betalist = []
    while beta_max - beta_min > deltabeta:
        beta = 0.5*(beta_min + beta_max)
        sis.beta = np.array([beta,beta])
        sis.randinit_upto([1,1])
        sis.x1 = np.power(sis.x1,0.2)
        sis.x2 = np.power(sis.x2,0.2)
        th1th2 = sis.integrate_by_chunks(tfinal,stop_condition=True)
        betalist.append(beta)
        outlist.append(list(np.vsplit(th1th2,2)[0].flatten()))
        final = outlist[-1][-1]
        print("%.7f -> %.8f"%(beta,final))
        if final > threshold:
            beta_max = beta
        else:
            beta_min = beta
    betac = beta_max
    for beta in [betac-np.logspace(1.1*np.log10(deltabeta),-2.5)]:
        sis.beta = np.array([beta, beta])
        sis.randinit_upto([1, 1])
        sis.x1 = np.power(sis.x1, 0.2)
        sis.x2 = np.power(sis.x2, 0.2)
        th1th2 = sis.integrate_by_chunks(tfinal, stop_condition=True)
        betalist.append(beta)
        outlist.append(list(np.vsplit(th1th2, 2)[0].flatten()))
        final = outlist[-1][-1]
        print("%.7f -> %.8f" % (beta, final))


json.dump([betalist,outlist],open("sisplateau_" + "_".join(argv[1:])+"_"+datetime.now().isoformat()+".json","w"))



