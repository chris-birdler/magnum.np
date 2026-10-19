import json, re, sys, numpy as np
sys.path.insert(0, ".")
from scan2 import main_dirs, RUNS
rows = []
for d in main_dirs() + sorted(RUNS.glob("TB_*_a5")):
    s = json.loads((d / "summary.json").read_text()); ini = s.get("init", {})
    h = ini.get("relax_history", [])
    if not h:
        continue
    E = np.array([x["E_over_Kd"] for x in h]); dE = np.diff(E)
    w9 = None
    st = [x for x in s["stages"] if x.get("measure")]
    if st:
        c = st[0]["cycles"][2:] or st[0]["cycles"][-1:]
        w9 = np.mean([x["w_loop"] for x in c])
    rows.append(dict(name=d.name, conv=ini.get("relax_converged"), it=ini.get("relax_iterations"), E=E[-1],
                     up_last2=int((dE[-2:] > 0).sum()), up_all=int((dE > 0).sum()), n=len(dE),
                     trms=h[-1]["torque_rms"], tmax=h[-1]["torque_max"], absdE=np.abs(dE[-2:]).max(), w9=w9,
                     ulp=np.abs(dE[-2:]).max() / (abs(E[-1]) * 2**-23)))
E = np.array([r["E"] for r in rows]); it = np.array([r["it"] for r in rows])
print("runs", len(rows), "converged", sum(r["conv"] for r in rows), "iterations median/max", np.median(it), it.max())
print("E/Kd final: min %.4g median %.4g max %.4g; |E| min %.3g; n E<0: %d" % (E.min(), np.median(E), E.max(), np.abs(E).min(), (E < 0).sum()))
print("runs with energy RISE in one of the 2 quiet chunks: %d; chunks with rise over all history: %d of %d" % (
    sum(r["up_last2"] > 0 for r in rows), sum(r["up_all"] for r in rows), sum(r["n"] for r in rows)))
u = np.array([r["ulp"] for r in rows]); print("max |dE| of the 2 quiet chunks in float32 units of E: median %.1f 90%% %.1f max %.1f" % (np.median(u), np.quantile(u, .9), u.max()))
tr = np.array([r["trms"] for r in rows]); tm = np.array([r["tmax"] for r in rows])
print("final torque rms A/m median %.2f max %.2f; torque max median %.1f max %.1f" % (np.median(tr), tr.max(), np.median(tm), tm.max()))
# allowed rest drift per cycle (6.67 chunks of 5 ns per 33.3 ns cycle) relative to the 9 mT loop loss
r9 = np.array([1e-6 * abs(r["E"]) * 6.67 / r["w9"] for r in rows if r["w9"]])
act = np.array([r["absdE"] * 6.67 / r["w9"] for r in rows if r["w9"]])
print("criterion bound: 1e-6|E| per chunk x 6.67 chunks/cycle / w_loop(9 mT): median %.3f 90%% %.3f max %.3f" % (np.median(r9), np.quantile(r9, .9), r9.max()))
print("actual last |dE| per chunk x 6.67 / w_loop(9 mT): median %.3f 90%% %.3f max %.3f" % (np.median(act), np.quantile(act, .9), act.max()))
for r in rows:
    if r["name"].startswith("TB"): print(r)
