from __future__ import division, print_function
import lambdac_results_loader as lrl
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cbook import flatten


def num(x):
    try:
        return int(x)
    except ValueError:
        return float(x)


k = num(argv[1]) if len(argv) > 1 else 12
dr = 0.17 if k == 50 else 0.37
"""
Plots the distance between the stable branch and the unstable branch in terms of r when starting from the
desynchronized state(1) the synchronized state(2) and d = lambda_f - lambda_b.
Requires db built in default location.
"""

l = lrl.Loader()
cur = l.con.cursor()
Nlist = cur.execute("select N from lamc where k=%i and t_final > 200 group by N" % k).fetchall()
Nlist = list(flatten(Nlist))
maxN = max(Nlist)

lv, fv, m = l.plot_jump_gaps(transpose=True, k=k)
f, lavgf, lerrf = l.get_lambdac_of_f(maxN, "lambda_f", k=k, where=" t_final > 200 ")
plt.plot(f, lavgf, lw=2, label="$\lambda_f$", color=l.color_list[0])
CS = plt.contour(fv, lv, m.T, levels=[dr], linewidths=2, color="k")
CS.collections[0].set_label("∆R = "+str(dr))
if k == 50:
    plt.text(0.5, 0.01, "desynchronized", horizontalalignment='center', verticalalignment='center')
    plt.text(0.5, 0.15, "synchronized", horizontalalignment='center', verticalalignment='center')
    plt.ylim(ymax=0.201)
elif k == 12:
    plt.text(0.5, 0.05, "desynchronized", horizontalalignment='center', verticalalignment='center')
    plt.text(0.4, 0.4, "synchronized", horizontalalignment='center', verticalalignment='center')
plt.gcf().set_size_inches(9.07,6.67,True)
plt.tight_layout()
plt.legend(loc="best")
plt.show()
plt.savefig("/tmp/lambda_forward_k%i.png"%k,dpi=600)


lv, fv, m = l.plot_jump_gaps_backwards(transpose=True, k=k)
f, lavgb, lerrb = l.get_lambdac_of_f(maxN, "lambda_b", k=k, where=" t_final > 200 ")
plt.plot(f, lavgb, lw=2, label="$\lambda_b$", color=l.color_list[0])
CS = plt.contour(fv, lv, m.T, levels=[dr], linewidths=2, color="k")
CS.collections[0].set_label("∆R = "+str(dr))
if k == 50:
    # plt.clabel(CS,inline=1,fontsize=16)
    plt.text(0.5, 0.01, "desynchronized", horizontalalignment='center', verticalalignment='center')
    plt.text(0.5, 0.15, "synchronized", horizontalalignment='center', verticalalignment='center')
    plt.ylim(ymax=0.201)
elif k == 12:
    CS = plt.contour(fv, lv, m.T, levels=[dr], linewidths=2, color="k")
    # plt.clabel(CS,inline=False,fontsize=16)
    plt.text(0.5, 0.05, "desynchronized", horizontalalignment='center', verticalalignment='center')
    plt.text(0.4, 0.4, "synchronized", horizontalalignment='center', verticalalignment='center')
plt.gcf().set_size_inches(9.07,6.67,True)
plt.tight_layout()
plt.legend(loc="best")
plt.show()

plt.savefig("/tmp/lambda_backward_k%i.png"%k,dpi=600)
l.plot_all_d(errorbar=True, k=k, where=" t_final > 190 ")
f,lfv,lbv = l.get_predicted_d(delta_r=dr,k=k)
plt.plot(f,lfv-lbv,lw=3,color="k",label="theory")
plt.ylim(ymin=0)
plt.gcf().set_size_inches(9.85,7.85,True)
plt.tight_layout()
plt.legend(loc="best")
plt.show()
try:
    plt.savefig("/tmp/f_of_d_k%i.pdf"%k)
except RuntimeError:
    plt.savefig("/tmp/f_of_d_k%i.png" % k,dpi=600)