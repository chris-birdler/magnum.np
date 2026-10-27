#!/usr/bin/env python
"""
Precision test of PLAN 13. For each patch, the arms old and new ran the same drive cycle with --precision_cycles double;
step_dump.py --trace saved m after 1, 2, 4, ..., 256 sample intervals, and after each sample interval the counts of
RKF45 tries and accepted steps, t and m in 2048 probe cells.
- ρ: the RMS over the magnetic cells (|m0| above 0) of |m_new − m_old|.
- Yardstick: ρ_null, the median ρ of the 3 S3 null runs in single precision (s3_eval.py JSON: null[i].div_rms_dm), the
  scatter of production runs that round-off alone causes.
- Criterion: ρ after the cycle in double precision at most BOUND · ρ_null (BOUND 1e-3; B1, model change: proposal 0.1).
- The same sequence of tries and accepted steps in both arms after each sample interval. If it differs first after
  sample interval j, the result is "not determined" and the probe ρ up to j − 1 is given.
- Positive control (--control): the old code with the drive amplitude changed by the relative amount REL. Detection
  limit: the relative drive change that gives ρ = BOUND · ρ_null (linear response).
- A pair with bitwise equal S3 end states in single precision passes as in S2 (--bitwise NAME).

    python prec_eval.py --pair NAME OLD NEW S3_JSON INIT BOUND [--pair ...] [--control NAME BASE CTRL REL S3_JSON INIT]
                        [--bitwise NAME ...] --json OUT
"""
import argparse
import json
import math
import statistics

import torch


def cycle(d):
    c = json.load(open(d + "/summary.json"))["stages"][0]["cycles"][0]
    return c["w_loop"], c["steps"]


def compare(old, new, init):
    A, B = torch.load(old + "/m.pt", map_location="cpu"), torch.load(new + "/m.pt", map_location="cpu")
    mag = torch.load(init, map_location="cpu", weights_only=False)["m"].double().norm(dim=-1) > 0
    course = []
    for k in sorted(A["m"]):
        ma, mb = A["m"][k].double(), B["m"][k].double()
        ang = torch.atan2(torch.linalg.cross(ma, mb).norm(dim=-1), (ma * mb).sum(-1))[mag]  # exact also near 0
        course.append({"k": k, "rho": float((mb - ma).norm(dim=-1)[mag].pow(2).mean().sqrt()),
                       "max_angle_deg": math.degrees(float(ang.max())), "dtype": str(A["m"][k].dtype)})
    ta, tb = A["trace"], B["trace"]
    assert torch.equal(A["probe_idx"], B["probe_idx"]), "different probe cells"
    n = min(len(ta["tries"]), len(tb["tries"]))
    j = next((i for i in range(n) if (ta["tries"][i], ta["accepted"][i]) != (tb["tries"][i], tb["accepted"][i])), None)
    probe = [float((pb.double() - pa.double()).norm(dim=-1).pow(2).mean().sqrt())
             for pa, pb in zip(ta["probe_m"], tb["probe_m"])]
    (w_old, s_old), (w_new, s_new) = cycle(old), cycle(new)
    return {"course": course, "probe_rho": probe, "first_step_difference": None if j is None else j + 1,
            "tries_old": ta["tries"][-1], "tries_new": tb["tries"][-1], "accepted_old": ta["accepted"][-1],
            "accepted_new": tb["accepted"][-1], "t_equal": ta["t"] == tb["t"], "steps_old": s_old, "steps_new": s_new,
            "D_new_double": math.log(w_new) - math.log(w_old), "rho_double": course[-1]["rho"],
            "float64": all(c["dtype"] == "torch.float64" for c in course)}


def null_rho(s3):
    r = json.load(open(s3))
    return statistics.median(x["div_rms_dm"] for x in r["null"]), r["new"]["div_rms_dm"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", nargs=6, action="append", default=[],
                    metavar=("NAME", "OLD", "NEW", "S3_JSON", "INIT", "BOUND"))
    ap.add_argument("--control", nargs=6, action="append", default=[],
                    metavar=("NAME", "BASE", "CTRL", "REL", "S3_JSON", "INIT"))
    ap.add_argument("--bitwise", nargs="*", default=[])
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    res = {"pairs": [], "controls": [], "bitwise_single": a.bitwise}
    for name, base, ctrl, rel, s3, init in a.control:
        r = compare(base, ctrl, init)
        rn, _ = null_rho(s3)
        r.update({"name": name, "rel_drive_change": float(rel), "rho_null_single": rn,
                  "detection_limit_rel_drive_1e-3": float(rel) * 1e-3 * rn / r["rho_double"] if r["rho_double"] > 0
                  else float("inf")})  # inf: the control changed nothing, the test does not see the change
        res["controls"].append(r)
        print("control %-7s drive × (1 + %.0e): ρ double %.2e; ρ_null single %.2e; detection limit for ρ ≤ 1e-3 ρ_null: "
              "relative drive change %.1e; first step difference %s; D double %+.2e"
              % (name, float(rel), r["rho_double"], rn, r["detection_limit_rel_drive_1e-3"],
                 r["first_step_difference"], r["D_new_double"]))
    for name, old, new, s3, init, bound in a.pair:
        bound = float(bound)
        r = compare(old, new, init)
        rn, rp = null_rho(s3)
        ratio = r["rho_double"] / rn
        if r["first_step_difference"] is not None:
            verdict = "not determined (tries differ first after sample interval %d)" % r["first_step_difference"]
        else:
            verdict = "pass" if ratio <= bound else "FAIL"
        r.update({"name": name, "rho_null_single": rn, "rho_pair_single": rp, "ratio": ratio, "bound": bound,
                  "result": verdict})
        res["pairs"].append(r)
        print("%-8s ρ double %.2e; ρ_null single %.2e (pair in single %.2e); ratio %.1e (bound %.0e): %s"
              % (name, r["rho_double"], rn, rp, ratio, bound, verdict))
        print("         tries old/new %d/%d, accepted %d/%d, t equal %s, steps in the cycle %d/%d, D double %+.2e%s"
              % (r["tries_old"], r["tries_new"], r["accepted_old"], r["accepted_new"], r["t_equal"], r["steps_old"],
                 r["steps_new"], r["D_new_double"], "" if r["float64"] else ", WARNING: m is not float64"))
        print("         ρ (largest angle in degrees) after sample interval " + ", ".join(
            "%d: %.1e (%.1e)" % (c["k"], c["rho"], c["max_angle_deg"]) for c in r["course"]))
        print("         probe ρ (2048 cells) after sample interval 1, 16, 64, 128, 256: " + ", ".join(
            "%.1e" % r["probe_rho"][i - 1] for i in (1, 16, 64, 128, 256) if i <= len(r["probe_rho"])))
    for name in a.bitwise:
        print("%-8s S3 end states in single precision bitwise equal: pass (as S2), no run in double precision" % name)
    json.dump(res, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
