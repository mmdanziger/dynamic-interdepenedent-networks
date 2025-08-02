from __future__ import division,print_function
import sqlite3
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import json
from glob import glob


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


def str_or_quote(x):
    if type(x) is str:
        return "\"" + x + "\""
    else:
        return str(x)


def add_individual_solution_to_db(intfunc, lam, k, db_path="/home/micha/phd/cudamoto_lambdac.db"):
    import single_net_sync as sns
    fracvec = [0.01 * i for i in range(101)]
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    for frac in fracvec:
        rsols = sns.solve_R(intfunc, lam, frac, k)
        elevels = [sns.get_energy_level(intfunc, lam, frac, k, i) for i in rsols]
        insert_statement = ""
        if len(rsols) == 1:
            insert_statement = "insert into sols (f,k,lambda) values(%.4f,%.4f,%.9f)" % (frac, k, lam)
        elif len(rsols) == 2:
            insert_statement = "insert into sols (f,k,lambda,rsol1,energy1) values(%.4f,%.4f,%.9f,%.9f,%.9f)" % (
                frac, k, lam, rsols[1], elevels[1])
        elif len(rsols) == 3:
            insert_statement = "insert into sols (f,k,lambda,rsol1,energy1,rsol2,energy2) values(%.4f,%.4f,%.9f,%.9f,%.9f,%.9f,%.9f)" % (
                frac, k, lam, rsols[1], elevels[1], rsols[2], elevels[2])
        print(insert_statement)
        cur.execute(insert_statement)
    cur.close()
    con.commit()
    con.close()


def add_lambda_c_runs_to_db(lambda_c_glob, dry=True, db_path="/home/micha/phd/cudamoto_lambdac.db"):
    con = sqlite3.connect(db_path)
    # con = sqlite3.connect(":memory:") #for debug purposes
    cur = con.cursor()
    skipped, added = [0, 0]
    for fname in glob(lambda_c_glob):
        k, v = zip(*json.load(open(fname)).items())
        wherestring = " and ".join([str(k_) + " = " + str_or_quote(v_) for k_, v_ in zip(k, v)])
        if cur.execute("SELECT count(*) FROM lamc WHERE " + wherestring).fetchone()[0]:
            skipped += 1
            continue
        added += 1
        if not dry:
            cur.execute("INSERT INTO lamc (" + ",".join(map(str, k)) + ") VALUES (" + ",".join(map(str, v)) + ")")
    print("Added %i, skipped %i" % (added, skipped))
    cur.close()
    con.commit()
    con.close()


def build_db(db_path="/home/micha/phd/cudamoto_lambdac.db", lambda_c_glob="/home/micha/lambdac/*.json",
             trhist_glob="/home/micha/phd/oscillators/processed/*",
             energy_fname="/home/micha/phd/oscillators/zhang_results_energy_HD.json"):
    con = sqlite3.connect(db_path)
    # con = sqlite3.connect(":memory:") #for debug purposes
    cur = con.cursor()
    cur.execute(
        "CREATE TABLE " +
        "lamc(lamid INTEGER PRIMARY KEY AUTOINCREMENT, N INT, f num, h num, k num, lambda_b num, lambda_f num, " +
        "precision num, r_threshold num, t_final INT, timestamp INT)")
    for fname in glob(lambda_c_glob):
        k, v = zip(*json.load(open(fname)).items())
        cur.execute("INSERT INTO lamc (" + ",".join(map(str, k)) + ") VALUES (" + ",".join(map(str, v)) + ")")
    cur.execute("CREATE TABLE trhist_run(trid INTEGER PRIMARY KEY AUTOINCREMENT, N INT, f num, k num, fname TEXT)")
    cur.execute(
        "CREATE TABLE " +
        "trhist_values(trvalueid INTEGER PRIMARY KEY AUTOINCREMENT, trid INTEGER, synchronized INT, lambda num, " +
        "r num, rstd num, rlineartrend num)")

    for tr_fname in glob(trhist_glob):
        d = json.load(open(tr_fname))
        d["N"] = 2 ** d["N"]
        tr_run_keys = ["N", "f", "k", "fname"]
        cur.execute("INSERT INTO trhist_run (" + ",".join(tr_run_keys) + ") VALUES(" + ",".join(
            [str_or_quote(d[k]) for k in tr_run_keys]) + ")")
        cur.execute("SELECT max(trid) FROM trhist_run")
        current_trid = cur.fetchone()[0]
        for row in d["synchronized"]:
            cur.execute("INSERT INTO trhist_values (trid,synchronized,lambda,r,rstd,rlineartrend) VALUES(" + ",".join(
                map(str_or_quote,
                    [current_trid, 1, row["lambda"], row["r"], row["rstd"], row["rtrends"][0][0]])) + ")")
        for row in d["unsynchronized"]:
            cur.execute("INSERT INTO trhist_values (trid,synchronized,lambda,r,rstd,rlineartrend) VALUES(" + ",".join(
                map(str_or_quote,
                    [current_trid, 0, row["lambda"], row["r"], row["rstd"], row["rtrends"][0][0]])) + ")")

    cur.execute(
        "CREATE TABLE " +
        "sols (solid INTEGER PRIMARY KEY AUTOINCREMENT, f num, k num, lambda num, " +
        "rsol1 num DEFAULT NULL, rsol2 DEFAULT NULL, energy1 DEFAULT NULL, energy2 DEFAULT NULL)")
    try:
        d2 = json.load(open(energy_fname))
    except TypeError:
        con.commit()
        con.close()
        return con

    for f_str in d2["r"].keys():
        for lidx, l in enumerate(d2["lambda_vec"]):
            if len(d2["r"][f_str][lidx]) == 1:
                cur.execute("insert into sols (f,k,lambda) values(%.4f,%.4f,%.9f)" % (float(f_str), d2["k"], l))
            elif len(d2["r"][f_str][lidx]) == 2:
                cur.execute("insert into sols (f,k,lambda,rsol1,energy1) values(%.4f,%.4f,%.9f,%.9f,%.9f)" % (
                    float(f_str), d2["k"], l, d2["r"][f_str][lidx][1], d2["energy_level"][f_str][lidx][1]))
            elif len(d2["r"][f_str][lidx]) == 3:
                cur.execute(
                    "insert into sols (f,k,lambda,rsol1,energy1,rsol2,energy2) " +
                    "values(%.4f,%.4f,%.9f,%.9f,%.9f,%.9f,%.9f)" % (
                        float(f_str), d2["k"], l, d2["r"][f_str][lidx][1], d2["energy_level"][f_str][lidx][1],
                        d2["r"][f_str][lidx][2],
                        d2["energy_level"][f_str][lidx][2]))

    con.commit()
    con.close()

    return con


def add_sols_to_db(db_path="/home/micha/phd/cudamoto_lambdac.db",
                   energy_fname="/home/micha/phd/oscillators/zhang_results_energy_k50_HD.json", k=50):
    con = sqlite3.connect(db_path)
    # con = sqlite3.connect(":memory:") #for debug purposes
    cur = con.cursor()
    d2 = json.load(open(energy_fname))
    if "k" not in d2.keys():
        if k:
            d2["k"] = k
    for f_str in d2["r"].keys():
        for lidx, l in enumerate(d2["lambda_vec"]):
            if len(d2["r"][f_str][lidx]) == 1:
                cur.execute("insert into sols (f,k,lambda) values(%.4f,%.4f,%.9f)" % (float(f_str), d2["k"], l))
            elif len(d2["r"][f_str][lidx]) == 2:
                cur.execute("insert into sols (f,k,lambda,rsol1,energy1) values(%.4f,%.4f,%.9f,%.9f,%.9f)" % (
                    float(f_str), d2["k"], l, d2["r"][f_str][lidx][1], d2["energy_level"][f_str][lidx][1]))
            elif len(d2["r"][f_str][lidx]) == 3:
                cur.execute(
                    "insert into sols (f,k,lambda,rsol1,energy1,rsol2,energy2) " +
                    "values(%.4f,%.4f,%.9f,%.9f,%.9f,%.9f,%.9f)" % (
                        float(f_str), d2["k"], l, d2["r"][f_str][lidx][1], d2["energy_level"][f_str][lidx][1],
                        d2["r"][f_str][lidx][2],
                        d2["energy_level"][f_str][lidx][2]))

    con.commit()
    con.close()

    return con


class Loader(object):
    def __init__(self, fname="/home/micha/phd/cudamoto_lambdac.db"):
        self.fname = fname
        self.con = sqlite3.connect(fname)
        self.Nlist = self.make_N_list()
        self.color_list = ["darkorange", "darksage", "purple", "dodgerblue", "darkslategray", "red", "mediumblue",
                           "green", "sienna"]

    def __del__(self):
        self.con.close()

    def make_N_list(self):
        cur = self.con.cursor()
        cur.execute("SELECT N FROM lamc GROUP BY N")
        toreturn = [i[0] for i in cur.fetchall()]
        cur.close()
        return toreturn

    def get_lambdac_of_f(self, N, which_lambda="lambda_f", k=12, where=' r_threshold < 0.1 '):
        cur = self.con.cursor()
        if where:
            qstring = "select f,%s from lamc where N = %i and k = %i and %s order by f" % (which_lambda, N, k, where)
            print(qstring + ";")
            cur.execute(qstring)
        else:
            cur.execute("select f,%s from lamc where N = %i and k = %i order by f" % (which_lambda, N, k))
        res = defaultdict(list)
        for f_, lamf_ in cur.fetchall():
            res["%.5f" % f_].append(lamf_)

        lam_f_avg = [np.mean(lam_f_vec[1]) for lam_f_vec in sorted(res.items(), key=lambda x: float(x[0]))]
        lam_f_std = [np.std(lam_f_vec[1]) for lam_f_vec in sorted(res.items(), key=lambda x: float(x[0]))]
        f = [float(lam_f_vec[0]) for lam_f_vec in sorted(res.items(), key=lambda x: float(x[0]))]
        cur.close()
        return [f, lam_f_avg, lam_f_std]

    def plot_r_of_lambda(self, k=12, flist=[0, 0.1, 0.3, 0.5, 0.7, 1]):
        cur = self.con.cursor()
        this_clist = self.color_list[:]
        for f in flist:
            cur.execute(
                "select lambda,rsol1,rsol2,energy1,energy2 from sols where abs(f - %.3f) < 0.001 and abs(k - %.3f) < 0.001 order by lambda asc" % (
                f, k))
            lamv_all, lamv_r1, lamv_r2, r1v, r2v = [[], [], [], [], []]
            founde2 = False
            result_set = list(cur.fetchall())
            for lam, r1, r2, e1, e2 in result_set:
                lamv_all.append(lam)

                if r1 and r2:
                    r1, r2 = sorted([r1, r2])
                    lamv_r1.append(lam)
                    r1v.append(r1)
                    lamv_r2.append(lam)
                    r2v.append(r2)
                elif r1 and not r2:
                    lamv_r2.append(lam)
                    r2v.append(r1)
                elif r2 and not r1:
                    lamv_r2.append(lam)
                    r2v.append(r2)
                if not founde2 and e2 and e2 < 0:
                    lam_center = lam
            c = this_clist.pop()
            plt.plot(lamv_all, [0 for i in lamv_all], color=c)
            plt.plot(lamv_r1, r1v, ":", lw=3, color=c)
            plt.plot(lamv_r2, r2v, "-", color=c, lw=3, label="f = %.2f" % f)
            # plt.axvline(lam_center,ls=":",color=c)
        cur.close()
        plt.xlabel("$\lambda$")
        plt.ylabel("R")
        plt.legend(loc="best")
        plt.show()

    def plot_all_d(self, errorbar=False, k=12, where=" t_final > 180 "):
        plt.figure()
        from matplotlib.colors import hsv_to_rgb
        hue_of_N = lambda N:  0.75*(np.log2(N) - np.log2(min(self.Nlist))) / (np.log2(max(self.Nlist)) -np.log2(min(self.Nlist)) )

        for N in reversed(self.Nlist):
            try:
                f, lavg, lerr = self.get_lambdac_of_f(N, "lambda_f - lambda_b", k=k, where=where)
            except:
                continue
            if not f:
                continue
            this_color = hsv_to_rgb((hue_of_N(N),1,1))
            if errorbar:
                plt.errorbar(f, lavg, color=this_color,ecolor=this_color, yerr=lerr, label="N = %i" % N)
            else:
                plt.plot(f, lavg, '.-', color=this_color,label="N = %i" % N)
            print("%i gets "%N,end=" ")
            print(this_color)
        plt.legend(loc="best")
        plt.xlabel("$f$")
        plt.ylabel("$d$")

    def plot_all_lambda(self, lambda_type="lambda_f", errorbar=False):
        for N in self.Nlist:
            f, lavg, lerr = self.get_lambdac_of_f(N, lambda_type)
            if errorbar:
                plt.errorbar(f, lavg, yerr=lerr, label="$N = %i$" % N)
            else:
                plt.plot(f, lavg, '.-', label="$N = %i$" % N)
        plt.legend(loc="best")
        plt.xlabel("$f$")
        plt.ylabel("$\\" + lambda_type + "$")

    def get_delta_r_at_lambda_c(self, lambda_type="lambda_f"):
        cur = self.con.cursor()
        cur.execute("select f,N,sum(%s)/count(*) from lamc group by f,N " % lambda_type)
        result_set = list(cur.fetchall())
        rstdvec, rsol_e_vec = [[], []]
        for f, N, lamc in result_set:
            cur.execute(
                "select rstd from trhist_values  join trhist_run using (trid) " +
                "where synchronized = 0 and abs( f - %.3f) < 1e-4 and N = %i order by abs(lambda - %.9f) asc limit 0,1" % (
                    f, N, lamc))
            res = cur.fetchone()
            if res:
                rstdvec.append(res[0])
            else:
                print("Missing rstd result for f,N,lambdac = " + ",".join(map(str, [f, N, lamc])))
                rstdvec.append(None)
            cur.execute(
                "select rsol1,energy1 from sols where abs( f - %.3f) < 1e-4 order by abs(lambda - %.9f) asc limit 0,1" % (
                    f, lamc))
            res = cur.fetchone()
            if res:
                rsol_e_vec.append(res)
            else:
                print("Missing sol for result for f,lambdac = " + ",".join(map(str, [f, lamc])))
                rsol_e_vec.append(None)
        joined = filter(lambda x: all(x), zip(result_set, rstdvec, rsol_e_vec))
        cur.close()
        return list(joined)

    def get_predicted_d(self,delta_r=0.37,k=50):
        cur = self.con.cursor()
        lambda_f_qstring = "select f,min(lambda) from sols where k=%i and rsol1>0 and rsol1 < %.9f group by f order by f;"%(k,delta_r)
        lambda_b_qstring = "select f,max(lambda) from sols where k=%i and rsol1>0 and rsol2 is not null and  (rsol2 - rsol1) < %.9f group by f order by f;"%(k,delta_r)
        cur.execute(lambda_f_qstring)
        res={}
        for f,lf in cur.fetchall():
            res["%.3f"%f] = [lf]
        cur.execute(lambda_b_qstring)
        for f,lb in cur.fetchall():
            try:
                res["%.3f"%f].append(lb)
            except KeyError:
                print("No lf found for %.3f"%f)
        minimum_lambda_b_string = "select f,min(lambda) from sols where k=%i and rsol1>0 and f>0 group by f;"%k
        cur.execute(minimum_lambda_b_string)
        for f,lbmin in cur.fetchall():
            if "%.3f"%f in res and len(res["%.3f"%f]) < 2:
                res["%.3f"%f].append(lbmin)

        fvec, llv = zip(*sorted(map(lambda x: (float(x[0]),x[1]), res.items())))
        lfvec,lbvec = zip(*llv)
        return np.array(fvec), np.array(lfvec), np.array(lbvec)

    def plot_jump_gaps(self, data_type="rsol1", k=12, transpose=False):
        from matplotlib import cm as colormaps
        cur = self.con.cursor()
        cur.execute("SELECT f FROM sols where abs(k - %.5f) <1e-4 GROUP BY f order by f" % k)
        fvec = np.array([i[0] for i in cur.fetchall()])
        cur.execute("SELECT lambda FROM sols where abs(k - %.5f) <1e-4  GROUP BY lambda order by lambda" % k)
        lamvec = np.array([i[0] for i in cur.fetchall()])
        cur.execute(
            "SELECT f,lambda,rsol1,energy1 FROM sols where abs(k - %.5f) <1e-4 group by lambda,f ORDER BY f,lambda " % k)
        A = np.zeros((len(fvec) * len(lamvec)))
        res = cur.fetchall()
        try:
            assert (len(res) == len(fvec) * len(lamvec))
        except AssertionError:
            print("lambda = %i , fvec = %i , res = %i" % (len(lamvec), len(fvec), len(res)))
            raise ValueError
        for idx, [f, lam, rgap, egap] in enumerate(res):
            if data_type == "rsol1":
                A[idx] = rgap if egap and egap > 0 else np.nan
            else:
                A[idx] = egap if egap and egap > 0 else np.nan
        A = A.reshape((len(fvec), len(lamvec)))
        mask = np.ma.masked_array(data=A, mask=np.isnan(A))
        plt.figure()
        if transpose:
            plt.pcolormesh(fvec, lamvec, mask.T, cmap=colormaps.viridis_r)
            plt.xlabel("$f$")
            plt.ylabel("$\lambda$")
        else:
            plt.pcolormesh(lamvec, fvec, mask, cmap=colormaps.viridis_r)
            plt.xlabel("$\lambda$")
            plt.ylabel("$f$")
        plt.colorbar()
        plt.show()
        return lamvec, fvec, mask

    def plot_jump_gaps_backwards(self, data_type="rsol1", k=12, transpose=True):
        from matplotlib import cm as colormaps
        cur = self.con.cursor()
        cur.execute("SELECT f FROM sols where abs(k - %.5f) <1e-4 GROUP BY f order by f" % k)
        fvec = np.array([i[0] for i in cur.fetchall()])
        cur.execute("SELECT lambda FROM sols where abs(k - %.5f) <1e-4  GROUP BY lambda order by lambda" % k)
        lamvec = np.array([i[0] for i in cur.fetchall()])
        cur.execute(
            "SELECT f,lambda,rsol1,rsol2,energy1,energy2 FROM sols where abs(k - %.5f) <1e-4 group by lambda,f ORDER BY f,lambda" % k)
        A = np.zeros((len(fvec) * len(lamvec)))
        for idx, [f, lam, r1, r2, e1, e2] in enumerate(cur.fetchall()):
            if data_type == "rsol1":
                A[idx] = r2 - r1 if e2 and e1 and e2 < e1 else np.nan
            else:
                A[idx] = e1 - e2 if e2 and e1 and e2 < e1 else np.nan
        A = A.reshape((len(fvec), len(lamvec)))
        mask = np.ma.masked_array(data=A, mask=np.isnan(A))
        plt.figure()
        if transpose:
            plt.pcolormesh(fvec, lamvec, mask.T, cmap=colormaps.viridis_r)
            plt.xlabel("$f$")
            plt.ylabel("$\lambda$")
        else:
            plt.pcolormesh(lamvec, fvec, mask, cmap=colormaps.viridis_r)
            plt.xlabel("$\lambda$")
            plt.ylabel("$f$")
        plt.colorbar()
        plt.show()
        return lamvec, fvec, mask

    def plot_relative_gap(self, k=12, transpose=True):
        from matplotlib import cm as colormaps
        cur = self.con.cursor()
        cur.execute("SELECT f FROM sols where abs(k - %.5f) <1e-4 GROUP BY f" % k)
        fvec = np.array([i[0] for i in cur.fetchall()])
        cur.execute("SELECT lambda FROM sols where abs(k - %.5f) <1e-4  GROUP BY lambda" % k)
        lamvec = np.array([i[0] for i in cur.fetchall()])
        cur.execute(
            "SELECT f,lambda,rsol1,rsol2,energy1,energy2 FROM sols where abs(k - %.5f) <1e-4 ORDER BY f,lambda" % k)
        A = np.zeros((len(fvec) * len(lamvec)))
        for idx, [f, lam, r1, r2, e1, e2] in enumerate(cur.fetchall()):
            A[idx] = (r2 - r1) / r1 if e2 and e1 and e2 < e1 else np.nan
        A = A.reshape((len(fvec), len(lamvec)))
        mask = np.ma.masked_array(data=A, mask=np.isnan(A))
        plt.figure()
        if transpose:
            plt.pcolormesh(fvec, lamvec, mask.T, cmap=colormaps.viridis_r)
            plt.xlabel("$f$")
            plt.ylabel("$\lambda$")
        else:
            plt.pcolormesh(lamvec, fvec, mask, cmap=colormaps.viridis_r)
            plt.xlabel("$\lambda$")
            plt.ylabel("$f$")
        plt.colorbar()
        plt.show()
        return lamvec, fvec, mask
