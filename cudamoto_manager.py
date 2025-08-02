from __future__ import division
import json
from os import path
from glob import glob
import numpy as np
from matplotlib import pyplot as plt

import numpy


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y


def num(x):
    try:
        return int(x)
    except ValueError:
        return float(x)


def load_rhist(fname):
    d = json.load(open(fname))
    l_vec, r_vec = zip(*d)
    out = {}
    basename = path.basename(fname)

    # basename = basename.strip(".json")
    # basename = basename[:basename.find(".", basename.find(".") + 1)]

    for entry in basename.split("_")[2:-1]:
        out[entry[0]] = num(entry[1:])
    out["lambda"] = l_vec
    out["r"] = r_vec
    return out


def get_lambda_c(rhist):
    lvec = rhist["lambda"]
    rvec = rhist["r"]
    synchronized = False
    lambda_f, lambda_b = [0, 0]
    for idx, (l, r) in enumerate(zip(lvec, rvec)):
        if r > 0.1 and not synchronized:
            synchronized = True
            lambda_f = l
        if r < 0.1 and synchronized:
            lambda_b = l
            synchronized = False
    return (lambda_f, lambda_b)


def load_glob(glob_string):
    results = [[], []]
    for fname in glob(glob_string):
        rhist = load_rhist(fname)
        lambda_c = get_lambda_c(rhist)
        results[1].append(
                {"f": rhist["f"], "lambda_f": lambda_c[0], "lambda_b": lambda_c[1], "d": lambda_c[0] - lambda_c[1]})
        results[0].append(rhist)
    results[0].sort(key=lambda x: x["f"])
    results[1].sort(key=lambda x: x["f"])
    return results


def plot_result_list(results):
    for rhist in results:
        plt.plot(rhist["lambda"], rhist["r"], ".-", label="$f = %.2f$" % rhist["f"])
    plt.legend(ncol=2, loc="best")


def plot_critical_values(critical_results):
    plt.figure()
    critical_vectors = dict((k, list(map(lambda x: x[k], critical_results))) for k in critical_results[0].keys())
    plt.plot(critical_vectors["f"], critical_vectors["d"], ".-", label="$d$")
    plt.plot(critical_vectors["f"], critical_vectors["lambda_f"], ".-", label=r"$\lambda_f$")
    plt.plot(critical_vectors["f"], critical_vectors["lambda_b"], ".-", label=r"$\lambda_b$")
    plt.axhline(y=0, ls='-.')
    plt.legend(ncol=2, loc="best")


def load_trhist(fname, check_trend_to=3, r_dividing=0.2):
    d = json.load(open(fname))
    out = {}
    basename = path.basename(fname)
    from scipy import polyfit
    # basename = basename.strip(".json")
    # basename = basename[:basename.find(".", basename.find(".") + 1)]
    for entry in basename.split("_")[2:-1]:
        out[entry[0]] = num(entry[1:])
    out["fname"] = basename
    out["synchronized"] = []
    out["unsynchronized"] = []
    for lambda_run in d:
        trends = []
        t, r = zip(*lambda_run["dynamics"])
        l = lambda_run["lambda"]
        rfinal = lambda_run["r"]
        if rfinal > r_dividing:
            run_key = "synchronized"
        else:
            run_key = "unsynchronized"
        rstd = np.std(r)
        for i in range(check_trend_to):
            trends.append(polyfit(t, r, i + 1).tolist())
        out[run_key].append({"lambda": l, "r": rfinal, "rstd": rstd, "rtrends": trends})
    return out


def load_trhist_glob(glob_string):
    results = []
    flist = glob(glob_string)
    for idx, fname in enumerate(flist):
        print("Loading file %i of %i (%s)" % (idx, len(flist), fname))
        results.append(load_trhist(fname))
    return results


def reduce_trhist(fname, out_dir="."):
    basename = path.basename(fname)
    ofname = path.join(out_dir, "processed" + basename)
    if path.isfile(ofname):
        return

    out = load_trhist(fname)
    json.dump(out, open(ofname, "w"))


def construct_dict_from_dir(globstring, trend_filter=0.001, r_cutoff=0.2):
    if trend_filter is None:
        return [json.load(open(fname)) for fname in glob(globstring)]
    else:
        output = []
        for fname in glob(globstring):
            this_run = json.load(open(fname))
            syncd = [x for x in this_run["synchronized"] if
                     abs(x["rtrends"][0][0]) < trend_filter and x["r"] > r_cutoff]
            unsyncd = [x for x in this_run["unsynchronized"] if
                       abs(x["rtrends"][0][0]) < trend_filter and x["r"] < r_cutoff]
            this_run["synchronized"] = syncd
            this_run["unsynchronized"] = unsyncd
            output.append(this_run)
    return output


def plot_all_of_size(results_list, N):
    results_list = list(filter(lambda x: x["N"] == N, results_list))
    f1 = plt.figure()
    ax1 = f1.add_subplot(111)
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)

    for result in sorted(results_list, key=lambda x: x["f"]):
        l, rstd = zip(*map(lambda x: (x["lambda"], x["rstd"]), result["synchronized"]))
        ax1.plot(l, rstd, label="$f = %.2f$" % result["f"])
        l, rstd = zip(*map(lambda x: (x["lambda"], x["rstd"]), result["unsynchronized"]))
        ax2.plot(l, rstd, label="$f = %.2f$" % result["f"])
    ax1.set_title("synchronized state")
    ax2.set_title("unsynchronized state")
    ax1.set_xlabel("$\lambda$")
    ax1.set_ylabel("$\sigma_r$")
    ax2.set_xlabel("$\lambda$")
    ax2.set_ylabel("$\sigma_r$")
    ax1.legend(loc="best")
    ax2.legend(loc="best")
    plt.show()


def plot_all_of_f(results_list, f, size_scaling_function=None):
    if size_scaling_function is None:
        size_scaling_function = lambda x: 1
    results_list = list(filter(lambda x: "%.3f" % x["f"] == "%.3f" % f, results_list))
    f1 = plt.figure()
    ax1 = f1.add_subplot(111)
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)

    for result in sorted(results_list, key=lambda x: x["N"]):
        l, rstd = zip(*map(lambda x: (x["lambda"], x["rstd"]), result["synchronized"]))
        ax1.plot(l, np.array(rstd) * size_scaling_function(result["N"]), label="$N = 2^{%i}$" % result["N"])
        l, rstd = zip(*map(lambda x: (x["lambda"], x["rstd"]), result["unsynchronized"]))
        ax2.plot(l, np.array(rstd) * size_scaling_function(result["N"]), label="$N = 2^{%i}$" % result["N"])
    ax1.set_title("synchronized state")
    ax2.set_title("unsynchronized state")
    ax1.set_xlabel("$\lambda$")
    ax1.set_ylabel("$\sigma_r$")
    ax2.set_xlabel("$\lambda$")
    ax2.set_ylabel("$\sigma_r$")
    ax1.legend(loc="best")
    ax2.legend(loc="best")
    plt.show()


def smooth_l_r(l, r, w):
    lr = sorted(zip(l, r))
    output = []
    for idx, lrpair in enumerate(lr):
        if idx > w and idx < len(lr) - w:
            avg_val = np.mean([x[1] for x in lr[idx - w:idx + w]])
            output.append((lrpair[0], avg_val))
    l, r = zip(*output)
    return [np.array(l), np.array(r)]


def phase_flipping_updown_count(run):
    t, r = zip(*run["dynamics"])
    current_state = 1
    states = []
    thresh = 0.4
    state_count = 0
    for val in r:
        if val > thresh:
            if val < 1.2 * thresh:
                continue
            if current_state == 1:
                state_count += 1
            else:
                states.append([current_state, state_count])
                current_state = 1
                state_count = 1
        else:
            if val > 0.8 * thresh:
                continue
            if current_state == 0:
                state_count += 1
            else:
                states.append([current_state, state_count])
                current_state = 0
                state_count = 0
    return states
