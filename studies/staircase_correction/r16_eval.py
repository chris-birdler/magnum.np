#!/usr/bin/env python
"""
PLAN 16: the release rules of a speed patch.

--t2 NAME NEW ABP: R2 and R3 (single precision; T2 files of s1_states.py --against, each tree against its own
double-precision output, 16 states) of a tree of class R against the released tree ABP.
  - R2 (accuracy): for each quantity the RMS over the states of e; deterministic quantities: e_new at most 1.10 × e_ABP;
    noise quantities (err, err_rms, norm): at most 1.25 × e_ABP. A field whose reference has an RMS below 1 cu in all
    states: the absolute RMS difference (key abs) at most 1e-3 cu more than ABP (rule of 14a).
  - R3 (systematic part, deterministic quantities): u = 2^-24, d_k = b_new,k - b_ABP,k, SE = sample SD / sqrt(n):
    pass if (a) |mean d| <= max(u, 4 SE(d)) (no coherent change above the noise and above u), or (b) |mean b_new| <=
    max(u, |mean b_ABP|) (the systematic part does not grow). For each quantity the script also gives the detection
    limit: the smallest shift of b in all states (both signs) that fails R3.
  Exit code 0 if all quantities pass, else 1.

--same A B: class E, bitwise: two T2 files (e, b, abs, ref and the repetitions of every quantity and state) or two T1
  folders of s1_states.py (every tensor and number of every state file state_*.pt, tensors compared as bytes, so +0.0
  and -0.0 differ; meta.pt holds the output paths and is not compared). Exit code 0 if all are equal, else 1.

--selftest T2FILE: the rules on a real T2 file and on changed copies: the file against itself passes; e × 1.2 of one
  deterministic quantity fails R2; a shift of b by 2u away from 0 in all states of one field fails R3; the same shift
  in one state only passes R3; a shift toward 0 that leaves |mean b| below u passes R3 by (b) (the systematic part
  falls); a shift through 0 to more than |mean b_ABP| on the other side fails R3.

    python r16_eval.py --t2 AS5 t2/B/AS5.pt t2/B/ABP.pt --json res/r16_B.json
    python r16_eval.py --same t1/B/AS3 t1/B/AS5
"""
import argparse
import copy
import json
import math
import pathlib
import statistics
import sys

import torch

NOISE = ("err", "err_rms", "norm")
U = 2.0 ** -24
E_DET, E_NOISE = 1.10, 1.25


def detection(md, mn, la, lb):
    """The smallest shift delta of b in all states (in units of u, one for each sign) that fails R3: |md + delta| > la
    and |mn + delta| > lb (a constant shift does not change SE(d)). Scan with a relative step of 0.2 %."""
    res = []
    for sg in (1.0, -1.0):
        x, hit = 1e-3, None
        while x < 1e7:
            dl = sg * x * U
            if abs(md + dl) > la and abs(mn + dl) > lb:
                hit = sg * x
                break
            x *= 1.002
        res.append(hit if hit is not None else sg * math.inf)
    return res


def rules(new, old, name="", out=print):
    """R2 and R3 of the tree new against old (loaded T2 records). Returns (pass, result per quantity)."""
    cu = 1e-6 * new["meta"]["H_amp"]
    ks = sorted(new["states"])
    if sorted(old["states"]) != ks:
        raise SystemExit("r16_eval: the two files have other states")
    res, ok = {}, True
    for q in new["states"][ks[0]]["q"]:
        if q not in old["states"][ks[0]]["q"]:
            out("     %-34s not in the ABP file: information only" % q)
            continue
        qn = [new["states"][k]["q"][q] for k in ks]
        qo = [old["states"][k]["q"][q] for k in ks]
        rms = lambda v: math.sqrt(sum(x * x for x in v) / len(v))
        small = q.startswith("h") and max(x.get("ref", 1e300) for x in qn) < cu
        if small or "e" not in qn[0]:
            an, ao = rms([x["abs"] for x in qn]) / cu, rms([x["abs"] for x in qo]) / cu
            r = {"rule": "abs", "abs_cu_new": an, "abs_cu_ABP": ao, "pass": an <= ao + 1e-3}
        else:
            en, eo = rms([x["e"] for x in qn]), rms([x["e"] for x in qo])
            lim = E_NOISE if q in NOISE else E_DET
            r = {"rule": "noise" if q in NOISE else "det", "e_new": en, "e_ABP": eo, "e_ratio": en / eo if eo else
                 (0.0 if en == 0 else math.inf), "e_bound": lim}
            r["R2"] = en <= lim * eo
            if q not in NOISE:
                bn, bo = [x["b"] for x in qn], [x["b"] for x in qo]
                d = [x - y for x, y in zip(bn, bo)]
                n = len(d)
                se_d = statistics.stdev(d) / math.sqrt(n)
                md, mn, mo = statistics.mean(d), statistics.mean(bn), statistics.mean(bo)
                la, lb = max(U, 4 * se_d), max(U, abs(mo))
                a_ok, b_ok = abs(md) <= la, abs(mn) <= lb
                r.update({"mean_d": md, "lim_d": la, "mean_b_new": mn, "lim_b_new": lb, "mean_b_ABP": mo,
                          "R3a": a_ok, "R3b": b_ok, "R3": a_ok or b_ok, "detect_u": detection(md, mn, la, lb)})
            r["pass"] = r["R2"] and r.get("R3", True)
        res[q] = r
        ok &= r["pass"]
    out("R2/R3 %-8s against ABP: %s (%d states)" % (name, "PASS" if ok else "FAIL", len(ks)))
    for q, r in res.items():
        if r["rule"] == "abs":
            txt = "abs %.2e cu (ABP %.2e cu)" % (r["abs_cu_new"], r["abs_cu_ABP"])
        else:
            txt = "e ratio %.4f (bound %.2f)" % (r["e_ratio"], r["e_bound"])
            if "R3" in r:
                txt += "; mean d %+.2e (limit %.2e)%s; mean b new %+.2e, ABP %+.2e%s; fails from %s u" % (
                    r["mean_d"], r["lim_d"], "" if r["R3a"] else " *", r["mean_b_new"], r["mean_b_ABP"],
                    "" if r["R3b"] else " *", "/".join("%+.1f" % x for x in r["detect_u"]))
        out("     %-34s %s%s" % (q, txt, "" if r["pass"] else "  <- FAIL"))
    return ok, res


def same(a, b, path=""):
    """Bitwise equality of two nested records; returns a list of the paths that differ."""
    if torch.is_tensor(a) or torch.is_tensor(b):
        if not (torch.is_tensor(a) and torch.is_tensor(b)) or a.dtype != b.dtype or a.shape != b.shape:
            return [path]
        ac, bc = a.detach().cpu().contiguous(), b.detach().cpu().contiguous()
        return [] if ac.numel() == 0 or torch.equal(ac.view(torch.uint8), bc.view(torch.uint8)) else [path]
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            return [path + " (keys)"]
        return [p for k in a for p in same(a[k], b[k], "%s/%s" % (path, k))]
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b):
            return [path + " (length)"]
        return [p for i, (x, y) in enumerate(zip(a, b)) for p in same(x, y, "%s[%d]" % (path, i))]
    if isinstance(a, float) and isinstance(b, float):
        return [] if (math.isnan(a) and math.isnan(b)) or a.hex() == b.hex() else [path]
    return [] if a == b else [path]


def load(p):
    return torch.load(p, map_location="cpu", weights_only=False)


def same_main(a, b):
    pa, pb = pathlib.Path(a), pathlib.Path(b)
    if pa.is_dir():
        fa, fb = sorted(x.name for x in pa.glob("state_*.pt")), sorted(x.name for x in pb.glob("state_*.pt"))
        if fa != fb or not fa:
            print("same: the folders have other files (%d, %d)" % (len(fa), len(fb)))
            return False
        diff = [p for f in fa for p in same(load(pa / f), load(pb / f), f)]
        n = len(fa)
    else:
        ra, rb = load(pa), load(pb)
        diff = same(ra["states"], rb["states"], "states")
        n = len(ra["states"])
    print("same %s %s: %d %s, %d differences%s" % (a, b, n, "files" if pa.is_dir() else "states", len(diff),
                                                   (": " + ", ".join(diff[:10])) if diff else ""))
    return not diff


def selftest(f):
    r = load(f)
    quiet = lambda *x: None
    ok_self, _ = rules(r, r, "self", quiet)
    ks = sorted(r["states"])
    det = [q for q in r["states"][ks[0]]["q"] if q not in NOISE and "e" in r["states"][ks[0]]["q"][q]
           and q.startswith("h")]
    q = det[0]
    e12 = copy.deepcopy(r)
    for k in ks:
        e12["states"][k]["q"][q]["e"] *= 1.2
    ok_e, res_e = rules(e12, r, "e x1.2", quiet)
    mb = statistics.mean(r["states"][k]["q"][q]["b"] for k in ks)
    sg = 1.0 if mb >= 0 else -1.0
    b2 = copy.deepcopy(r)
    for k in ks:
        b2["states"][k]["q"][q]["b"] += sg * 2 * U
    ok_b, res_b = rules(b2, r, "b+2u", quiet)
    b1 = copy.deepcopy(r)
    b1["states"][ks[0]]["q"][q]["b"] += sg * 2 * U
    ok_b1, res_b1 = rules(b1, r, "b+2u one", quiet)
    b0 = copy.deepcopy(r)
    for k in ks:
        b0["states"][k]["q"][q]["b"] -= mb - 0.5 * sg * U
    ok_b0, res_b0 = rules(b0, r, "b to 0.5u", quiet)
    bx = copy.deepcopy(r)
    for k in ks:
        bx["states"][k]["q"][q]["b"] -= 2.5 * mb
    ok_bx, res_bx = rules(bx, r, "b through 0", quiet)
    checks = {"file against itself passes": ok_self,
              "e x 1.2 of %s fails R2" % q: (not ok_e) and not res_e[q]["R2"],
              "b shifted 2u away from 0 in all states of %s fails R3" % q: (not ok_b) and not res_b[q]["R3"],
              "the same shift in one state of %s passes R3" % q: res_b1[q]["R3"],
              "b shifted to mean 0.5u in %s passes R3 by (b)" % q: res_b0[q]["R3b"] and not res_b0[q]["R3a"],
              "b shifted through 0 to -1.5 x its mean in %s fails R3" % q: not res_bx[q]["R3"]}
    for k, v in checks.items():
        print("selftest: %-50s %s" % (k, "ok" if v else "WRONG"))
    return all(checks.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t2", nargs=3, metavar=("NAME", "NEW", "ABP"))
    ap.add_argument("--same", nargs=2, metavar=("A", "B"))
    ap.add_argument("--selftest", metavar="T2FILE")
    ap.add_argument("--json")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest(a.selftest) else 1)
    if a.same:
        sys.exit(0 if same_main(*a.same) else 1)
    if a.t2:
        ok, res = rules(load(a.t2[1]), load(a.t2[2]), a.t2[0])
        if a.json:
            pathlib.Path(a.json).write_text(json.dumps({"name": a.t2[0], "new": a.t2[1], "ABP": a.t2[2], "pass": ok,
                                                        "quantities": res}, indent=1))
        sys.exit(0 if ok else 1)
    ap.error("one of --t2, --same, --selftest")


if __name__ == "__main__":
    main()
