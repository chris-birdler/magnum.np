"""Tests of PLAN 4.21 (gap test of V1 with the production line): the box chain chain421.sh with a fake dispatcher
(order, set-up gate, first-run check, the stops R4 and D5), the parser of watch421.py, and the K5 end-to-end test
(the production argument set with V1 and E3 on a small geometry with contacts, then the checks and the
post-processing of the chain). The K5 test on the CPU is gate 1 of the plan (slow: set RUN_K5=1); on a GPU it is the
box test B5. Also the harmonic split of acc_reduce with the moment weight on the K5 geometry."""
import json
import math
import os
import pathlib
import shlex
import subprocess
import sys
import threading
import time

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "box421"))


# --- the chain with a fake dispatcher ---------------------------------------------------------------------------

FAKE_PY = r"""#!/bin/bash
# fake python of the chain test: check_sat_run.py prints PASS (FAIL if the run has R4FAIL_TEST); instance_postproc.py
# writes POSTPROC_DONE of the run of its --pattern; the other scripts do nothing
case "$1" in
  check_sat_run.py) if [ -f "$2/R4FAIL_TEST" ]; then echo "FAIL x"; echo "C5 data part: FAIL"; else echo "C5 data part: PASS"; fi ;;
  instance_postproc.py) n=$(echo "$@" | sed 's/.*--pattern \^\(.*\)\$.*/\1/'); touch "$CHAIN_B/runs/$n/POSTPROC_DONE" ;;
  *) exit 0 ;;
esac
"""
FAKE_SUBMIT = r"""#!/bin/bash
f=$CHAIN_B/spool/$1_$2.jobs
cp "$3" $f.tmp && mv $f.tmp $f && echo "submitted $f"
"""


class FakeDispatcher(threading.Thread):
    """One job per GPU (n_gpus), smaller priority first, then the name (as run_queue.py); a job writes stdout.log, the
    set-up end line after t_setup, then DONE (and a Bloch sample if the name is in bloch) after t_run; the spool file
    gets .done. STOP_NEW_JOBS: no new job starts."""

    def __init__(self, B, n_gpus=4, t_setup=1.0, t_run=None, bloch=(), r4fail=(), nodone=()):
        super().__init__(daemon=True)
        self.B, self.n, self.t_setup, self.t_run = pathlib.Path(B), n_gpus, t_setup, t_run or {}
        self.bloch, self.r4fail, self.nodone = set(bloch), set(r4fail), set(nodone)
        self.started, self.stop = [], False

    def run_job(self, f, name):
        d = self.B / "runs" / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "stdout.log").write_text("[tie] {}\n")
        time.sleep(self.t_setup)
        with open(d / "stdout.log", "a") as fh:
            fh.write("[relax] virgin state: converged=True after 0 iterations, 1 s wall\n")
        rows = "# cycle,sample,n_bloch_pos,n_bloch_neg\n0,0,0,0\n0,1,%d,0\n" % (1 if name in self.bloch else 0)
        (d / "bloch_00_amp00.csv").write_text(rows)        # during the run, as the online count of run_loops
        time.sleep(self.t_run.get(name, 2.0))
        if name in self.r4fail:
            (d / "R4FAIL_TEST").write_text("x")
        if name not in self.nodone:                 # a failed job ends without DONE
            (d / "DONE").write_text("ok\n")
        f.with_name(f.name + ".done").write_text(json.dumps({"done": [name], "failed": {}}))

    def run(self):
        spool, busy = self.B / "spool", []
        while not self.stop:
            busy = [t for t in busy if t.is_alive()]
            if len(busy) < self.n and not (self.B / "runs" / "STOP_NEW_JOBS").exists():
                for f in sorted(spool.glob("*.jobs"), key=lambda p: (int(p.name.split("_")[0]), p.name)):
                    name = f.read_text().split()[0]
                    if name in self.started or f.with_name(f.name + ".done").exists():
                        continue
                    self.started.append(name)
                    t = threading.Thread(target=self.run_job, args=(f, name), daemon=True)
                    t.start()
                    busy.append(t)
                    break
            time.sleep(0.2)


def _chain(tmp_path, d5="a", **kw):
    B = tmp_path / "B"
    (B / "jobs").mkdir(parents=True)
    (B / "spool").mkdir()
    (B / "runs").mkdir()
    (B / "jobs" / "plan421.txt").write_text((HERE / "jobs" / "plan421.txt").read_text())
    py, sub = tmp_path / "fakepy.sh", tmp_path / "submit.sh"
    py.write_text(FAKE_PY)
    sub.write_text(FAKE_SUBMIT)
    py.chmod(0o755)
    sub.chmod(0o755)
    disp = FakeDispatcher(B, **kw)
    disp.start()
    env = dict(os.environ, CHAIN_B=str(B), CHAIN_SUBMIT=str(sub), CHAIN_POLL="0.3", CHAIN_PY=str(py))
    r = subprocess.run(["bash", str(HERE / "chain421.sh"), "setup", d5], env=env, cwd=str(B), timeout=300,
                       capture_output=True, text=True)
    disp.stop = True
    log = (B / "runs" / "chain421.log").read_text()
    return r, disp, log, (B / "runs" / "CHAIN421_END").read_text(), B


def test_chain421_order_and_gates(tmp_path):
    """All 12 jobs run; the small runs start first; no large run starts before the first R4 PASS; each job is submitted
    only after the set-up of the one before ended; each finished run gets R4 PASS and the post-processing."""
    import conv_jobs as CJ
    t_run = {n + "_v16": (3.0 if n[1] in "13" else 4.0) for n in CJ.PLAN421_ORDER}
    r, disp, log, end, B = _chain(tmp_path, t_run=t_run)
    assert r.returncode == 0, r.stderr
    assert disp.started == [n + "_v16" for n in CJ.PLAN421_ORDER]
    assert "stopped" not in end and "12 jobs ended" in end
    lines = log.splitlines()
    first_pass = min(i for i, x in enumerate(lines) if "R4: C5 data part: PASS" in x)
    first_large = min(i for i, x in enumerate(lines) if "submitted K2_0_v16" in x)
    assert first_pass < first_large
    for n in disp.started:
        assert "%s: set-up ended" % n in log
        assert (B / "runs" / n / "POSTPROC_DONE").exists() and (B / "runs" / n / "R4_RESULT").read_text() == "PASS\n"


def test_chain421_stops(tmp_path):
    """D5 (a): a Bloch sample in a run of K1 stops new jobs; the jobs in progress end. An R4 FAIL stops new jobs, and
    the run keeps its files (no post-processing)."""
    r, disp, log, end, B = _chain(tmp_path / "d5", bloch={"K1_0_v16"})
    assert r.returncode == 0 and "D5 (a): a Bloch sample in K1_0_v16" in end
    assert len(disp.started) < 12 and (B / "runs" / "STOP_NEW_JOBS").exists()
    r, disp, log, end, B = _chain(tmp_path / "d5b", d5="b", bloch={"K1_0_v16"})
    assert r.returncode == 0 and "stopped" not in end and len(disp.started) == 12
    r, disp, log, end, B = _chain(tmp_path / "r4", r4fail={"K1_0_v16"})
    assert r.returncode == 0 and "R4 FAIL in K1_0_v16" in end
    assert not (B / "runs" / "K1_0_v16" / "POSTPROC_DONE").exists()
    assert "K2_0_v16" not in disp.started


def test_chain421_stop_with_a_queued_job(tmp_path):
    """A stop while a submitted job waits in the queue (2 GPUs): the job never starts, and the chain still ends. If
    all small runs end without DONE, no run passes R4: the chain stops new jobs and ends."""
    r, disp, log, end, B = _chain(tmp_path / "q", n_gpus=2, bloch={"K1_0_v16"}, t_run={"K1_0_v16": 4.0})
    assert r.returncode == 0 and "D5 (a)" in end
    import re
    waiting = re.findall(r"(K\d_\d_v16): not started before the stop", log)
    assert waiting and not set(waiting) & set(disp.started)
    small = {"K1_0_v16", "K3_0_v16", "K1_1_v16", "K3_1_v16"}
    r, disp, log, end, B = _chain(tmp_path / "n", nodone=small)
    assert r.returncode == 0 and "no small run passed R4" in end
    assert set(disp.started) == small


def test_watch421_parse():
    """watch421.parse_state: Bloch samples per run and file, the third pass, the GPU and host memory, the end."""
    import watch421 as W
    out = ("x\n2026-11-07 R4 FAIL in K1_0_v16\n@@STATUS\n{\"running\": {}}\n@@BLOCH\n"
           "/r/runs/K1_0_v16/bloch_00_amp00.csv 0\n/r/runs/K4_1_v16/bloch_06_amp06.csv 12\n@@PASS\n"
           "/r/runs/K2_0_v16/stdout.log\n@@MEM\n0, 15000\n1, 29000\ncg 100000000000 60000000000\nalive 1 0\n"
           "@@END\nchain421 end\n@@DONE\n3\n")
    s = W.parse_state(out)
    assert s["bloch"] == {"K4_1_v16": {"bloch_06_amp06.csv": 12}}
    assert s["pass3"] == ["K2_0_v16"] and s["gpu_mib"] == {0: 15000, 1: 29000} and s["host_avail_gb"] == 40.0
    assert s["alive"] == (1, 0)
    assert W.parse_state(out.replace("cg 100000000000", "cg max"))["host_avail_gb"] is None
    assert s["end"] == "chain421 end" and s["n_done"] == 3 and "R4 FAIL" in s["log"]


# --- K5: the production argument set on a small geometry with contacts -----------------------------------------

def k5_line():
    """The K4 line of the production (sat_corner_lines) at d/l_ex 48, phi 0.55, L_eff/l_ex 6, r_p 1.5 (gap about 1.6
    cells: the cut cells of two particles touch), 2 amplitudes (9 and 50 mT, snapshots at both), at most 3 cycles,
    plus the arguments of plan421_lines (the correction at q 16, V1, E3)."""
    import conv_jobs as CJ
    t = shlex.split(dict(CJ.sat_corner_lines(1))["K4_0_sat"])
    for a, v in (("d_lex", ["48.0"]), ("phi", ["0.55"]), ("Leff_lex", ["6.0"]), ("r_p", ["1.5"]), ("n_amp", ["2"]),
                 ("b_list", ["0.006", "0.033333"]), ("cycles_per_amp", ["2"]), ("max_cycles_per_amp", ["3"]),
                 ("snap_mT", ["9.0", "50.0"])):
        t = CJ.set_arg(t, a, v)
    return " ".join(t) + " " + CJ.PLAN421_ARGS


def check_k5_folder(d):
    """The readings of gate 1 / B5 on a finished K5 folder (after run_loops): R4 FAIL lines only for the stage count;
    the post-processing steps of chain421.sh end without an error; the results are consistent."""
    import domains
    cfg = json.loads((d / "config.json").read_text())
    assert cfg["n_touch"] > 0 and cfg["tie_cells"] == "thin_face" and cfg["diag"]["tie"]["n_tied"] > 0
    g = domains.run_geometry(d)
    g2 = domains.geometry_from_config(cfg)
    for k in ("ids", "frac", "inside"):
        assert np.array_equal(g[k], g2[k]), k
    assert np.allclose(g["mu"], g2["mu"], rtol=1e-6, atol=1e-9)
    assert list(d.glob("bloch_pos_*.csv")) and list(d.glob("diag_*.csv"))
    r = subprocess.run([sys.executable, str(HERE / "check_sat_run.py"), str(d)], capture_output=True, text=True)
    fails = [x for x in r.stdout.splitlines() if x.startswith("FAIL")]
    allowed = ("FAIL 7 stages complete", "FAIL 7 stage states", "FAIL snapshot stages")   # the test has 2 stages
    assert all(x.startswith(allowed) for x in fails), fails
    for script, extra in (("eval_hotspots.py", ["--json", str(d / "hotspots.json")]),
                          ("eval_hotcores.py", ["--json", str(d / "hotcores.json")]),
                          ("census_bloch.py", ["--out", str(d / "census_bloch.json")])):
        r = subprocess.run([sys.executable, str(HERE / script)] + extra + [str(d)], capture_output=True, text=True)
        assert r.returncode == 0, (script, r.stderr[-2000:])
    r = subprocess.run([sys.executable, str(HERE / "instance_postproc.py"), "--runs", str(d.parent), "--stop",
                        str(d.parent / "no_stop"), "--once", "--workers", "1", "--pattern", "^%s$" % d.name],
                       capture_output=True, text=True)
    assert (d / "POSTPROC_DONE").exists(), r.stdout[-2000:] + r.stderr[-2000:]
    acc = json.loads((d / "acc.json").read_text())
    for st in acc["stages"]:
        assert st["check_ok"]
        h = st["harm"]
        assert abs(sum(h["R_k"]) + h["R_rest"] - 1.0) < 1e-9     # an identity (R_rest = 1 - sum_k R_k): arithmetic only
        # Bessel's inequality (sum_k p_k <= p in each cell) only as information: the measured cycles are not exactly
        # periodic, and the corner runs K<c>_<k>_sat (voxel, no tie) have bessel_max above 1.01 in 104 of 259 stages
        # (PROTOKOLL 7.83; gate 1 of 4.21, eval_gate421_v1.py). The moment weight of p_k: test_acc_reduce_weights_*
        print("%s %s: bessel_max %.4f, cells above 1.001 %.4f" % (d.name, st["stage"], h["bessel_max"],
                                                                  h["bessel_frac_over_1e-3"]))
        assert abs(sum(p["p_share"] for p in st["particles"]) - 1.0) < 1e-9
        assert abs(sum(st["where"][k]["vol_share"] for k in ("r<0.5", "0.5<=r<0.8", "r>=0.8")) - 1.0) < 1e-9
    hs = json.loads((d / "hotspots.json").read_text())
    for run in hs if isinstance(hs, list) else [hs]:
        for stage in run["stages"].values():
            for name in ("nonf", "f"):
                for reg in stage[name]["regions"]:
                    assert isinstance(reg["particle"], int)
    return acc


@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_k5_end_to_end(tmp_path, device):
    """Gate 1 (CPU, RUN_K5=1) and box test B5 (GPU 0): run_loops with the K5 line, then the checks and the
    post-processing of chain421.sh on its folder (check_k5_folder)."""
    import torch
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("no GPU")
    if device == "cpu" and not os.environ.get("RUN_K5"):
        pytest.skip("slow: set RUN_K5=1 (gate 1 of PLAN 4.21)")
    d = tmp_path / "K5_v16"
    env = dict(os.environ, CUDA_DEVICE="-1" if device == "cpu" else os.environ.get("CUDA_DEVICE", "0"))
    r = subprocess.run([sys.executable, str(HERE / "run_loops.py")] + shlex.split(k5_line()) + ["--out", str(d)],
                       env=env, cwd=str(HERE), capture_output=True, text=True, timeout=5400)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    check_k5_folder(d)


# --- acc_reduce: the moment weight of the harmonic losses ---------------------------------------------------------

def _k5_geometry_cfg():
    """The geometry entries of config.json of the K5 runs of gate 1 (24^3 cells, phi 0.55, the fraction surface with
    the correction at q 16, V1): 9856 particle cells, 1824 of them tied, as in those runs."""
    dx, N = 1e-8, 24
    return {"N_used": N, "d": N * dx * (6 * 0.55 / (4 * math.pi)) ** (1 / 3), "dx": dx, "dx_lex": 3.0,
            "surface": "fraction", "staircase": "pair", "staircase_q": 16, "tie_cells": "thin_face", "ms_void": 1e-8}


def _harm_folder(d, cfg, g, mu_p, mu_file):
    """A run folder for acc_reduce.run with one stage of an exactly periodic motion: in each cell a constant plus the
    harmonics k 1, 2, 3 and 7, with an amplitude factor of the cell (0.5 ... 1.5). run_loops.HarmAcc accumulates m and
    the local loss p = alpha mu0 M_s / gamma |dm/dt|^2 mu_p (run_loops: p_coef p_w |m x h|^2 with p_w = M_s mu, the same
    loss, test_dissipation_formulas_agree); the acc file has the arrays of run_loops (--acc --acc_harm 7); geometry.npz
    (run_loops.write_run_geometry) has the moment weight mu_file."""
    import torch
    import run_loops as RL
    import acc_reduce as AR
    K, n, N = 7, 256, cfg["N_used"]
    d.mkdir()
    full = dict(cfg, alpha=0.1, Js=1.5, direction=[1.0, 0.0, 0.0], Leff_lex=6.0)
    (d / "config.json").write_text(json.dumps(full))
    RL.write_run_geometry(d / "geometry.npz", g, mu_file, N, cfg["d"], cfg["dx"], 1, True)
    mask = g["ids"] >= 0
    w, Ms = 2 * math.pi * 3.0e7, full["Js"] / AR.MU0
    c = np.random.default_rng(5).uniform(0.5, 1.5, mask.shape)
    amps = [(1, 0, 0.30, 0.2), (2, 1, 0.04, 0.3), (3, 1, 0.05, -1.0), (7, 2, 0.01, 2.0)]
    h = RL.HarmAcc((N, N, N, 3), K, "cpu", torch.as_tensor(mask))
    for j in range(n):
        th = 2 * math.pi * (j + 1) / n
        m = np.zeros((N, N, N, 3)) + np.array([0.1, -0.2, 0.3])
        dm = np.zeros_like(m)
        for k, comp, a, ph in amps:
            m[..., comp] += c * a * math.cos(k * th + ph)
            dm[..., comp] -= c * a * k * math.sin(k * th + ph)
        p = full["alpha"] * AR.MU0 * Ms / AR.GAMMA * w * w * (dm ** 2).sum(-1) * mu_p
        h.add(torch.as_tensor(m), torch.as_tensor(p), th, math.cos(th), float(p.mean()))
    h.end_cycle(True, 2)
    mt, s = h.mask, h.s_n
    np.savez(d / "acc_00_amp00.npz", p=(h.s["p"][mt] / s).float().numpy(), m0=(h.s["m0"][mt] / s).float().numpy(),
             m1=((h.s["re"][mt] + 1j * h.s["im"][mt]) * (2.0 / s)).to(torch.complex64).numpy(),
             h1=np.complex128(h.s_h1 * 2.0 / s), n_samples=s, cycles_from=2,
             p_box=float(h.s["p"].sum()) / s / h.s["p"].numel(), p_dis_samples=h.s_psum / s, freq=w / (2 * math.pi),
             mk2=(h.s_mk2[:, mt] / len(h.kept)).T.contiguous().numpy(), cycles_kept=np.array(h.kept),
             cycles_acc=h.n_cycles)
    (d / "samples_00_amp00.csv").write_text("# t,...,B\n" + "".join(",".join(["0"] * 7 + [b, "0"]) + "\n"
                                                                  for b in ("0.009", "-0.009")))


def test_acc_reduce_weights_the_harmonics_as_p(tmp_path):
    """The harmonic split of acc_reduce on the fraction surface with V1 (PLAN 4.21 K6): the loss p of the run carries
    the moment weight mu (run_loops gives the same tensor to the dissipation and to geometry.npz), and acc_reduce
    multiplies the harmonic losses p_k with mu of geometry.npz. For an exactly periodic motion with harmonics up to K
    only, sum_k p_k = p in each cell (Parseval): in the free cut cells and in the masters as in the full cells. Thus the
    largest ratio sum_k p_k / p is 1 and the rest is 0, in the box and in each particle. Controls: a weight that is
    missing on one side gives a largest ratio above 1 by more than 100 times the tolerance."""
    pytest.importorskip("torch")
    import domains
    import acc_reduce as AR
    cfg = _k5_geometry_cfg()
    g = domains.geometry_from_config(cfg)
    mask, mu, frac = g["ids"] >= 0, g["mu"], g["frac"]
    assert int(mask.sum()) == 9856 and int((mask & (mu == 0)).sum()) == 1824
    assert (mask & (mu > frac * 1.001)).any() and (mask & (frac < 1) & (mu > 0) & (mu < 1)).any()  # masters, cut cells
    mu_no_tie = domains.geometry_from_config(dict(cfg, tie_cells=None))["mu"]
    tol = 1e-5

    def harm(name, mu_p, mu_file):
        d = tmp_path / name
        _harm_folder(d, cfg, g, mu_p, mu_file)
        AR.run(d)
        st = json.loads((d / "acc.json").read_text())["stages"][0]
        assert st["check_ok"]
        bmax, rest = st["harm"]["bessel_max"], [st["harm"]["R_rest"]] + [q["R_rest"] for q in st["particles"]]
        print("%s: bessel_max - 1 %.3e, |R_rest| box %.2e, largest of a particle %.2e" % (
            name, bmax - 1.0, abs(rest[0]), max(abs(x) for x in rest[1:])))
        return bmax, rest

    bmax, rest = harm("right", mu, mu)
    assert abs(bmax - 1.0) < tol and max(abs(x) for x in rest) < tol, (bmax, rest)
    bmax, _ = harm("evaluation_without_mu", mu, np.ones_like(mu))    # p_k as on the voxel surface: ratio 1 / mu
    assert bmax > 1.0 + 100 * tol, bmax
    bmax, _ = harm("run_without_tie", mu_no_tie, mu)                 # masters: ratio group moment / f
    assert bmax > 1.0 + 100 * tol, bmax


# --- eval_421: the rules on synthetic run folders -----------------------------------------------------------------

def _fake_run(root, name, bloch=None, passes=1, last_bp=0, calm=(True,) * 7, r4="PASS", n_stages=7, pos=None,
              summary=True):
    """A run folder for eval_421.run_reading: summary.json (passes, stages), bloch_<i>_amp<i>.csv (bloch: {stage index:
    [(cycle, sample, n)]}), bloch_pos (pos: {stage index: [(cycle, sample, i, j, k)]}), geometry.npz (8^3 cells, one
    particle: cells 1 ... 6 inside, the shell of cells 0 and 7 a cut layer), R4_RESULT."""
    d = root / name
    d.mkdir(parents=True)
    pp = [{"bloch_points": [[0, 0, 0, 1, 0]] if (k < passes - 1 or last_bp) else [], "wall_s": 3600.0}
          for k in range(passes)]
    st = [{"name": "amp%02d" % i, "complete": True, "calm_ok": calm[i], "cycles": [{"wall_s": 600.0}]}
          for i in range(n_stages)]
    if summary:                                    # without: a run that stopped in the preparation (no summary.json)
        (d / "summary.json").write_text(json.dumps({"init": {"predemag_passes": pp, "bloch_free": last_bp == 0,
                                                             "wall_s": 60.0}, "stages": st}))
    for i in range(n_stages):
        rows = ["# cycle,sample,n_bloch_pos,n_bloch_neg"] + ["0,%d,0,0" % s for s in range(10)]
        for c, s, n in (bloch or {}).get(i, []):
            rows[1 + s] = "%d,%d,%d,0" % (c, s, n)
        (d / ("bloch_%02d_amp%02d.csv" % (i, i))).write_text("\n".join(rows) + "\n")
        prow = ["# cycle,sample,i,j,k,Q,particle,max_angle_cube_deg,f_edge_lo,f_edge_hi,f_cube_min"]
        prow += ["%d,%d,%d,%d,%d,1,0,170,0.5,0.5,0.5" % x for x in (pos or {}).get(i, [])]
        (d / ("bloch_pos_%02d_amp%02d.csv" % (i, i))).write_text("\n".join(prow) + "\n")
    ids = np.zeros((8, 8, 8), np.int8)
    inside = np.zeros((8, 8, 8), bool)
    inside[1:7, 1:7, 1:7] = True
    np.savez(d / "geometry.npz", ids=ids, frac=np.ones((8, 8, 8), np.float32), mu=np.ones((8, 8, 8), np.float32),
             inside=inside)
    (d / "R4_RESULT").write_text(r4 + "\n")
    (d / "DONE").write_text("ok\n")
    return d


def _verdict(tmp_path, special):
    import eval_421 as E
    R = []
    for c in (1, 2, 3, 4):
        for k in range(3):
            n = "K%d_%d_v16" % (c, k)
            R.append(E.run_reading(_fake_run(tmp_path, n, **special.get(n, {}))))
    return E.verdict(R), R


def test_eval421_rules(tmp_path):
    """The verdict of PLAN 4.21 on synthetic folders: clean -> passes; a Bloch sample at K2 is information; at K1 ->
    fails; no clean state after 3 passes at K3 -> fails (at K2: information); an R4 FAIL or 3 stages without a value
    -> not determined. The cube classes: a cube with all 8 cells inside is interior, a cube with a cell of the shell is
    in the cut layer; the longest run of consecutive cut-layer samples."""
    (v, why, nv), _ = _verdict(tmp_path / "a", {})
    assert v == "the production line with V1 passes" and nv == 0
    (v, why, nv), R = _verdict(tmp_path / "b", {"K2_0_v16": {"bloch": {6: [(0, 3, 1), (0, 4, 1), (0, 5, 1)]},
                                                             "pos": {6: [(0, 3, 0, 0, 0), (0, 4, 0, 0, 0),
                                                                         (0, 5, 3, 3, 3)]}},
                                               "K2_1_v16": {"passes": 3, "last_bp": 1}})
    assert v == "the production line with V1 passes"
    st = [s for r in R if r["run"] == "K2_0_v16" for s in r["R1"] if s["stage"] == "amp06"][0]
    assert (st["bloch_samples"], st["interior_samples"], st["cut_layer_samples"], st["longest_cut_layer_run"]) == \
        (3, 1, 2, 2)
    (v, why, nv), _ = _verdict(tmp_path / "c", {"K1_2_v16": {"bloch": {0: [(0, 1, 2)]}, "pos": {0: [(0, 1, 2, 2, 2)]}}})
    assert v == "the production line with V1 fails" and "K1_2_v16" in why[0]
    (v, why, nv), _ = _verdict(tmp_path / "d", {"K3_0_v16": {"passes": 3, "last_bp": 1}})
    assert v == "the production line with V1 fails"
    (v, why, nv), _ = _verdict(tmp_path / "e", {"K4_1_v16": {"r4": "FAIL"}})
    assert v == "not determined"
    (v, why, nv), _ = _verdict(tmp_path / "f", {"K2_2_v16": {"calm": (True,) * 4 + (False,) * 3}})
    assert v == "not determined" and nv == 3
    (v, why, nv), _ = _verdict(tmp_path / "g", {"K1_0_v16": {"n_stages": 4}, "K1_1_v16": {"n_stages": 5}})
    assert v == "not determined"
    (v, why, nv), _ = _verdict(tmp_path / "h", {"K4_0_v16": {"summary": False, "n_stages": 0}})
    assert v == "not determined" and any("R2 not evaluable in K4_0_v16" in x for x in why)
    (v, why, nv), _ = _verdict(tmp_path / "i", {"K1_2_v16": {"passes": 2, "last_bp": 1, "n_stages": 0}})
    assert v == "not determined" and any("R2 not evaluable in K1_2_v16" in x for x in why)
