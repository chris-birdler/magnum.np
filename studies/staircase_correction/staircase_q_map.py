#!/usr/bin/env python
"""
Static staircase error of the stray-field energy with the correction at q (PROTOKOLL 3.33), in % of the anisotropy
energy, for the points of the main β powder study. Law and data as in the deck staircase_energy_map (make_deck.py):
E_st error ≈ C (Δx/d)^p K_d from the B0b state at cells 3, 1.5, 0.75 l_ex (Richardson), |E_K| ≈ A_K Q_eff K_d.
The correction gives the exact field of the q-staircase (N3, PROTOKOLL 3.26 … 3.29), so its error is taken as that of
cells Δx/q (an estimate: the law is extrapolated to cells below 0.75 l_ex; m stays on the coarse cells).

    cd ~/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak && python <this script>
"""
import math, os, json, re, sys
import numpy as np
STUDY = os.path.expanduser("~/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak")
sys.path.insert(0, STUDY)
from anisotropy_model import box_edge
def md_rows(name):
    out = []
    for line in open(os.path.join(STUDY, "results", name)):
        if line.startswith("| ") and not line.startswith("| seed") and not line.startswith("|---"):
            out.append([c.strip() for c in line.strip().strip("|").split("|")])
    return out
def cfg(run): return json.load(open(os.path.join(STUDY, "runs", run, "config.json")))
def terms(rows, state, dx): return np.array([[float(r[3]), float(r[4]), float(r[5])] for r in rows if r[1] == state and float(r[2]) == dx])
ST = md_rows("mesh_energy_states.md")
B0 = {0.75: terms(ST, "F", 0.75), 1.5: terms(ST, "F→dx", 1.5), 3.0: terms(ST, "F→dx", 3.0)}
D_B0, Q_B0 = cfg("B0b_F_s941")["d_lex_used"], cfg("B0b_F_s941")["anisotropy"]["Q_eff"]
DM = {dx: B0[dx][:, 2].mean() for dx in B0}
R_ = (DM[3.0] - DM[1.5]) / (DM[1.5] - DM[0.75]); P = math.log2(R_)
E_LIM = DM[0.75] - (DM[1.5] - DM[0.75]) / (R_ - 1)
C = (DM[0.75] - E_LIM) / (0.75 / D_B0) ** P
A_K = abs(B0[0.75][:, 1].mean()) / Q_B0
def job_args(f):
    out = {}
    for line in open(os.path.join(STUDY, "jobs", f)):
        if line.strip() and not line.startswith("#"):
            n, _, a = line.strip().partition(" ")
            out[n] = {k: float(v) for k, v in re.findall(r"--(d_lex|Leff_lex|r_p|phi) ([0-9.]+)", a)}
    return out
CORN = job_args("sobol_corners.txt"); MAIN = job_args("sobol_main.txt")
K_D = 1.5**2 / (2*4e-7*math.pi)
print("order p %.2f, Q_B0 %.4g, D_B0 %.1f" % (P, Q_B0, D_B0))
P_Q = -math.log2(0.4845)          # reduction with q from the increments of X3 (PROTOKOLL 3.27): first order
def pct(x, Q, q, pq=None): return 100 * C * x ** P / (A_K * Q) / q ** (P if pq is None else pq)
# Q_eff from L_eff: K_eff = A / L_eff^2 ... use the deck relation? use Q = (l_ex/L_eff)^2 / 2? check below
for n in ("K1", "K3"):
    v = CORN[n + "_0"]; d = box_edge(v["d_lex"], v["phi"], 12.0)[1]; x = 3.0 / d
    Q = (1.0 / v["Leff_lex"]) ** 2   # K_eff/K_d = (l_ex/L_eff)^2 (A/K_eff = L_eff^2, l_ex^2 = A/K_d)
    print(n, "d/l_ex %.0f  Δx/d %.4f  Q %.5f (K_eff %.2f kJ/m3)" % (d, x, Q, Q * K_D / 1e3),
          " ".join("q%d %.0f%%" % (q, pct(x, Q, q)) for q in (1, 2, 4, 8, 16, 32)),
          "| with (1+r_p): " + " ".join("q%d %.0f%%" % (q, pct(x, Q * (1 + v["r_p"]), q)) for q in (1, 2, 4, 8, 16)))
# main study range
vals = []
for n, v in MAIN.items():
    d = box_edge(v["d_lex"], v.get("phi", 0.65), 12.0)[1]; x = 3.0 / d; Q = (1.0 / v["Leff_lex"]) ** 2
    vals.append([pct(x, Q * (1 + v.get("r_p", 0.0)), q) for q in (1, 2, 4, 8, 16)])
vals = np.array(vals)
print("main study (%d points, with (1+r_p)): " % len(vals) + " ".join("q%d %.0f … %.0f%% (median %.0f)" % (q, vals[:, i].min(), vals[:, i].max(), np.median(vals[:, i])) for i, q in enumerate((1, 2, 4, 8, 16))))
print("share of main points below 10 %: " + " ".join("q%d %.0f%%" % (q, 100 * np.mean(vals[:, i] < 10)) for i, q in enumerate((1, 2, 4, 8, 16))))
d = box_edge(96.0, 0.55, 12.0)[1]; x = 3.0 / d; Q = 1.0 / 144
print("staircase test geometry d/l_ex %.1f Δx/d %.4f Q %.5f: " % (d, x, Q) + " ".join("q%d %.0f%%" % (q, pct(x, Q, q)) for q in (1, 2, 4, 8)))

print("\nwith the reduction q^-%.2f of the X3 increments (3.27) instead of q^-%.2f:" % (P_Q, P))
vals2 = []
for n_, v in MAIN.items():
    d = box_edge(v["d_lex"], v.get("phi", 0.65), 12.0)[1]; x = 3.0 / d; Q = (1.0 / v["Leff_lex"]) ** 2
    vals2.append([pct(x, Q * (1 + v.get("r_p", 0.0)), q, P_Q) for q in (1, 2, 4, 8, 16)])
vals2 = np.array(vals2)
print("main study: " + " ".join("q%d median %.0f%% (%.0f … %.0f)" % (q, np.median(vals2[:, i]), vals2[:, i].min(), vals2[:, i].max()) for i, q in enumerate((1, 2, 4, 8, 16))))
print("share of main points below 10 %: " + " ".join("q%d %.0f%%" % (q, 100 * np.mean(vals2[:, i] < 10)) for i, q in enumerate((1, 2, 4, 8, 16))))
for n_ in ("K1", "K3"):
    v = CORN[n_ + "_0"]; d = box_edge(v["d_lex"], v["phi"], 12.0)[1]; x = 3.0 / d; Q = (1.0 / v["Leff_lex"]) ** 2 * (1 + v["r_p"])
    print(n_, "with (1+r_p): " + " ".join("q%d %.0f%%" % (q, pct(x, Q, q, P_Q)) for q in (1, 2, 4, 8, 16)))
