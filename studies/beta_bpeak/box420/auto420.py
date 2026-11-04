#!/usr/bin/env python3
"""PLAN 4.20, automatic box (go of Chris 2026-11-04 on PLAN 4.20, option 2, cost limit 16 $): rent a V100 box with
watch_rent420.rent_one(), test it and destroy it at once when a box test B1 ... B4 fails; then rent the next one, until
the failed boxes cost 1 $. A fail of B5 (the E6 tests on GPU 0) is a code fault: logs, destroy, stop, report. A box that
passes all tests starts the dispatcher (no automatic rerun: --retries 0 --gate_maxiter 0) and chain420.sh (stage B, the selection, stage
C). Then this script exits.
Box tests (PLAN 4.20, gates 3, as B1 ... B5 of PLAN 4.19):
  B1 running and ssh within 30 min (direct port of the container first, then the vast proxy);
  B2 GPU load: no compute process, utilization at most 5 %, the GPU count of the offer;
     fp32 matrix product at least 11.5 TFLOPS and copy bandwidth at least 700 GB/s on each GPU (free V100: 13.5, 824);
  B3 upload of the start states (D15_s941 ... 944, HFP_s921, 923, 924, 925, K4_2_sat for the speed check) with equal
     md5;
  B4 set-up box420_setup.sh <commit> (clone, pip, worktree) and the speed check off156 on each GPU at most 10 %
     slower than 3.699 ms per stage (box 11), in 75 min; a GPU more than 10 % faster is logged and passes;
  B5 test 6 of test_plan420.py on GPU 0 (end-to-end V1 and V2, a bitwise repeat).
Usage: auto420.py <commit> <gate A variants, e.g. V1,V2>. Destroy only instances that this script rented. Log:
auto420.log; state (spent $, failed machines): auto420_state.json."""
import hashlib, json, pathlib, re, subprocess, sys, time

import watch_rent420 as W

HERE = pathlib.Path(__file__).resolve().parent
LOG = HERE / "auto420.log"
STATE = HERE / "auto420_state.json"
BUDGET = 1.0                    # $ for failed boxes (PLAN 4.20, gates 3)
CREDIT_MIN = 17.5               # $ before a rent: the cost limit 16 $ plus 1.5 $ (PLAN 4.20, gates 8)
D15 = "/home/christoph/Documents/code/magnum.np-staircase/studies/staircase_correction/remote_20261020_box1/runs"
HFP = ("/home/christoph/Documents/code/magnum.np-diag/studies/beta_bpeak/runs/plan418_step3_box/x/magnum.np/studies/"
       "beta_bpeak/runs")
# (local file, place on the box, md5); the md5 values from d15_init.md5 (staircase box 6), plan418_step3_md5.txt, 4.18
UPLOADS = [(D15 + "/D15_s941/init.pt", "runs/D15_s941/init.pt", "ef9ecaf85ad270feb2053cf3a71fde8f"),
           (D15 + "/D15_s942/init.pt", "runs/D15_s942/init.pt", "db8c68b8ce4d11e34cc1e1457ec604ce"),
           (D15 + "/D15_s943/init.pt", "runs/D15_s943/init.pt", "308785e15d1b661da6f5df0a63480cff"),
           (D15 + "/D15_s944/init.pt", "runs/D15_s944/init.pt", "365a96f74c656e9bc304be6152771c28")]
UPLOADS += [(HFP + "/HFP_s%d/init.pt" % s, "runs/HFP_s%d/init.pt" % s, None) for s in (921, 923, 924, 925)]
INIT = pathlib.Path("/home/christoph/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak/runs/K4_2_sat/init.pt")
INIT_MD5 = "45bbbfa74a7b34c64b1d77bd3e1c9ab6"
FILES = [HERE / "box420_setup.sh", HERE / "gpucheck420.py"]
STAGE = "/root/magnum.np/studies/beta_bpeak/"
SSHO = ["-o", "StrictHostKeyChecking=no", "-o", "UserKnownHostsFile=/dev/null", "-o", "ConnectTimeout=15",
        "-o", "BatchMode=yes", "-o", "ServerAliveInterval=30"]
TFLOPS_MIN, GBS_MIN = 11.5, 700.0


def log(msg):
    line = time.strftime("%Y-%m-%d %H:%M:%S UTC ", time.gmtime()) + msg
    with LOG.open("a") as f:
        f.write(line + "\n")
    print(line, flush=True)


def load_state():
    return json.loads(STATE.read_text()) if STATE.exists() else {"spent": 0.0, "failed_machines": [], "boxes": []}


def save_state(s):
    STATE.write_text(json.dumps(s, indent=1))


def inst(iid):
    try:
        r = subprocess.run([W.V, "show", "instance", str(iid), "--raw"], capture_output=True, text=True, timeout=90)
        return json.loads(r.stdout)
    except Exception:          # also subprocess.TimeoutExpired: the caller tries again later
        return None


def credit():
    """The credit of the account in $, or None on an error."""
    try:
        r = subprocess.run([W.V, "show", "user", "--raw"], capture_output=True, text=True, timeout=90)
        return float(json.loads(r.stdout)["credit"])
    except Exception:
        return None


def instance_ids():
    """ids of all instances of the account, or None on an error"""
    try:
        r = subprocess.run([W.V, "show", "instances-v1", "--raw"], capture_output=True, text=True, timeout=90)
        j = json.loads(r.stdout)
        return [int(x["id"]) for x in j["instances"]] if j.get("success") else None
    except Exception:
        return None


def targets(d):
    """ssh targets of the container: the direct port of 22/tcp (or direct_port_start), then the vast proxy.
    Never machine_dir_ssh_port: that is the sshd of the host, not of our container."""
    t = []
    ports = d.get("ports") or {}
    p22 = ports.get("22/tcp") if isinstance(ports, dict) else None
    if p22 and d.get("public_ipaddr"):
        t.append((d["public_ipaddr"].strip(), int(p22[0]["HostPort"])))
    elif d.get("direct_port_start") and d.get("public_ipaddr"):
        t.append((d["public_ipaddr"].strip(), int(d["direct_port_start"])))
    if d.get("ssh_host") and d.get("ssh_port"):
        t.append((d["ssh_host"], int(d["ssh_port"])))
    return t


def ssh(tg, cmd, timeout=60):
    try:
        r = subprocess.run(["ssh"] + SSHO + ["-p", str(tg[1]), "root@" + tg[0], cmd], capture_output=True,
                           text=True, timeout=timeout)
        return r.returncode, r.stdout
    except subprocess.TimeoutExpired:
        return 124, ""


def scp_to(tg, src, dst, timeout):
    try:
        return subprocess.run(["scp", "-q"] + SSHO + ["-P", str(tg[1]), str(src), "root@%s:%s" % (tg[0], dst)],
                              capture_output=True, text=True, timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        return 124


def scp_from(tg, src, dst, timeout):
    try:
        return subprocess.run(["scp", "-q"] + SSHO + ["-P", str(tg[1]), "root@%s:%s" % (tg[0], src), str(dst)],
                              capture_output=True, text=True, timeout=timeout).returncode
    except subprocess.TimeoutExpired:
        return 124


def t1_ssh(iid):
    """Running and ssh in 30 min. Returns the working target or None."""
    t0, started = time.time(), False
    while time.time() - t0 < 1800:
        d = inst(iid)
        if d:
            st = d.get("actual_status")
            if st in ("stopped", "exited") and not started and time.time() - t0 > 300:
                log("T1 instance %s is %s: vastai start" % (iid, st))
                subprocess.run([W.V, "start", "instance", str(iid)], capture_output=True, text=True, timeout=90)
                started = True
            if st == "running":
                for tg in targets(d):
                    rc, out = ssh(tg, "echo SSH_OK", timeout=40)
                    if "SSH_OK" in out:
                        log("T1 pass: ssh %s:%d after %.0f s" % (tg[0], tg[1], time.time() - t0))
                        return tg
        time.sleep(30)
    log("T1 FAIL: no ssh in 30 min")
    return None


def t2_gpu(tg, n_gpus):
    rc, out = ssh(tg, "nvidia-smi --query-gpu=index,utilization.gpu --format=csv,noheader,nounits; echo APPS; "
                      "nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader", timeout=60)
    if rc != 0 or "APPS" not in out:
        log("T2 FAIL: nvidia-smi rc %d" % rc)
        return False
    head, apps = out.split("APPS", 1)
    try:
        util = [int(x.split(",")[1]) for x in head.strip().splitlines() if "," in x]
    except ValueError:
        log("T2 FAIL: utilization not readable: %r" % head[:200])
        return False
    apps = [a for a in apps.strip().splitlines() if a.strip()]
    log("T2 GPUs %d (offer %d), utilization %s, compute processes %d" % (len(util), n_gpus, util, len(apps)))
    if len(util) != n_gpus or apps or max(util) > 5:
        log("T2 FAIL: GPU count, foreign process or load")
        return False
    for f in FILES:
        if scp_to(tg, f, "/root/", 120) != 0:
            log("T2 FAIL: upload of %s" % f.name)
            return False
    rc, out = ssh(tg, "python3 /root/gpucheck420.py", timeout=300)
    res = re.findall(r"gpu (\d+) .*?: ([\d.]+) TFLOPS fp32, (\d+) GB/s", out)
    for g, tf, bw in res:
        log("T2 gpu %s: %s TFLOPS, %s GB/s" % (g, tf, bw))
    if len(res) != n_gpus or any(float(tf) < TFLOPS_MIN or float(bw) < GBS_MIN for _, tf, bw in res):
        log("T2 FAIL: throughput (limits %.1f TFLOPS, %.0f GB/s)" % (TFLOPS_MIN, GBS_MIN))
        return False
    log("T2 pass")
    return True


def t3_upload(tg):
    """B3: K4_2_sat (speed check) to /root/states, the start states to /root/upload/<place> (moved into the runs folder
    after the clone of B4); each md5 equal to the local file and to the list."""
    ssh(tg, "mkdir -p /root/states/K4_2_sat " + " ".join("/root/upload/" + str(pathlib.Path(p).parent)
                                                          for _, p, _ in UPLOADS), timeout=40)
    t = time.time()
    ok = True
    for src, dst, md5 in [(str(INIT), "/root/states/K4_2_sat/init.pt", INIT_MD5)] + [
            (a, "/root/upload/" + b, c) for a, b, c in UPLOADS]:
        loc = hashlib.md5(pathlib.Path(src).read_bytes()).hexdigest()
        rc = scp_to(tg, src, dst, 600)
        _, out = ssh(tg, "md5sum %s" % dst, timeout=60)
        good = rc == 0 and out.split()[:1] == [loc] and (md5 is None or md5 == loc)
        log("T3 %s %s: scp rc %d, md5 %s (local %s)" % (dst, "ok" if good else "FAIL", rc, out.strip()[:32], loc))
        ok &= good
    log("T3 %s in %.0f s" % ("pass" if ok else "FAIL", time.time() - t))
    return ok


def t4_setup(tg, n_gpus):
    ssh(tg, "cd /root; setsid nohup bash /root/box420_setup.sh %s > /root/setup420.out 2>&1 < /dev/null & disown"
        % COMMIT, timeout=40)
    t0, out = time.time(), ""
    while time.time() - t0 < 4500:
        time.sleep(60)
        rc, o = ssh(tg, "cat /root/progress.txt", timeout=40)
        if rc == 0 and o:
            out = o
        if "DONE setup420" in out or re.search(r"failed rc|no /root/states", out):
            break
    speed = re.findall(r"SPEED gpu (\d+) .*?-> (PASS|FAIL)", out)
    for line in out.splitlines():
        if "SPEED" in line or "failed" in line:
            log("T4 " + line.strip())
    ok = "DONE setup420" in out and len(speed) == n_gpus and all(v == "PASS" for _, v in speed)
    log("T4 %s after %.0f s" % ("pass" if ok else "FAIL", time.time() - t0))
    return ok


def t5_tests(tg):
    """B5: test 6 of test_plan420.py on GPU 0 (in 60 min)."""
    rc, out = ssh(tg, "cd %s; mkdir -p /root/b5; CUDA_DEVICE=-1 timeout 3600 python -m pytest -q -p no:logging "
                      "'test_plan420.py::test_run_loops_tie_end_to_end[0]' --basetemp /root/b5/tmp "
                      "> /root/b5/b5.log 2>&1; echo RC $?; tail -3 /root/b5/b5.log" % STAGE, timeout=3900)
    for line in out.strip().splitlines():
        log("T5 " + line)
    ok = "RC 0" in out and re.search(r"\b1 passed", out) is not None
    log("T5 %s" % ("pass" if ok else "FAIL"))
    return ok


def start_chain(tg, gate, mode):
    """Move the start states into runs/, start the GPU log, the dispatcher without retries and chain420.sh."""
    cmds = ("cd /root; cp -r /root/upload/runs/. %sruns/ && ls %sruns/; "
            "sed -i 's/--retries 1/--retries 0 --gate_maxiter 0/' /root/dispatcher.sh; "
            "grep -c 'retries 0 --gate_maxiter 0' /root/dispatcher.sh; "
            "setsid nohup nvidia-smi --query-gpu=timestamp,index,clocks.sm,temperature.gpu,power.draw,"
            "utilization.gpu,clocks_event_reasons.active --format=csv -l 60 >> /root/gpu_log.csv 2>&1 < /dev/null "
            "& disown; setsid nohup bash /root/dispatcher.sh >> /root/disp.out 2>&1 < /dev/null & disown; sleep 20; "
            "cd %s; setsid nohup bash chain420.sh %s %s > /root/chain420.out 2>&1 < /dev/null & disown; sleep 90; "
            "tail -5 %sdispatcher.log; cat %sruns/chain420.log" % (STAGE, STAGE, STAGE, gate, mode, STAGE, STAGE))
    rc, out = ssh(tg, cmds, timeout=240)
    for line in out.strip().splitlines():
        log("start " + line)
    return "submitted" in out


def fail_box(st, iid, d0, tg, why, offer, t_rent):
    """Small logs (md5), destroy, check, add the cost (from the instance data, else from the offer and the rent time)."""
    dst = HERE / ("box%s" % iid)
    dst.mkdir(exist_ok=True)
    if tg:
        ssh(tg, "cd /root; tar czf /root/diag%s.tgz progress.txt setup420.out speed416/*.json speed416/*.log "
                "speed416/*.csv b5/b5.log 2>/dev/null; md5sum /root/diag%s.tgz" % (iid, iid), timeout=60)
        if scp_from(tg, "/root/diag%s.tgz" % iid, dst, 120) == 0:
            f = dst / ("diag%s.tgz" % iid)
            log("logs %s md5 %s" % (f, hashlib.md5(f.read_bytes()).hexdigest()))
    gone = False
    for k in range(5):         # destroy until the instance is gone (5 tries); else the caller stops
        try:
            r = subprocess.run([W.V, "destroy", "instance", str(iid), "-y"], capture_output=True, text=True,
                               timeout=120)
            log("destroy %s (try %d) rc %d: %s" % (iid, k + 1, r.returncode, r.stdout.strip()[:120]))
        except subprocess.TimeoutExpired:
            log("destroy %s (try %d): timeout" % (iid, k + 1))
        time.sleep(15)
        ids = instance_ids()
        gone = ids is not None and int(iid) not in ids
        if gone:
            break
    d = dict(d0 or {})
    for k in ("machine_id", "host_id", "dph_total"):
        d[k] = d.get(k) or offer.get(k)
    hours = (time.time() - min(float(d.get("start_date") or t_rent), t_rent)) / 3600
    cost = hours * float(d["dph_total"])
    st["spent"] += cost
    if d.get("machine_id"):
        st["failed_machines"].append(int(d["machine_id"]))
    st["boxes"].append({"id": iid, "machine": d.get("machine_id"), "host": d.get("host_id"), "why": why,
                        "hours": round(hours, 3), "cost": round(cost, 3), "destroyed_check": gone})
    save_state(st)
    log("FAILED box %s (machine %s, %s): %.2f h, %.3f $; spent %.3f of %.2f $; instance gone %s"
        % (iid, d.get("machine_id"), why, hours, cost, st["spent"], BUDGET, gone))
    return gone


def main():
    global COMMIT
    COMMIT, gate = sys.argv[1], sys.argv[2]
    st = load_state()
    log("auto420 start: commit %s, gate A variants %s, budget %.2f $ for failed boxes, spent %.3f, failed machines %s"
        % (COMMIT, gate, BUDGET, st["spent"], st["failed_machines"]))
    while st["spent"] < BUDGET:
        cr = credit()
        if cr is None or cr < CREDIT_MIN:
            log("STOP: credit %s below %.2f $ (PLAN 4.20 gates 8: the cost limit plus 1.5 $)" % (cr, CREDIT_MIN))
            print("RESULT NO_CREDIT %s" % cr)
            return 5
        r = W.rent_one(excluded=st["failed_machines"])
        if r is None:
            log("STOP: no offer")
            print("RESULT NO_OFFER")
            return 1
        iid, offer = r
        t_rent = time.time()
        ram_gb = offer.get("cpu_ram", 0) * 2**20 / 1e9
        mode = "setup" if ram_gb >= 128.0 else "serial"
        log("rented %s: machine %s host %s, %d GPUs, %.3f $/h, host memory %.0f GB (q16 mode %s)" % (
            iid, offer.get("machine_id"), offer.get("host_id"), offer["num_gpus"], offer["dph_total"], ram_gb, mode))
        time.sleep(20)
        d0, tg, why = inst(iid), None, None
        try:
            tg = t1_ssh(iid)
            d0 = inst(iid) or d0
            if tg is None:
                why = "B1 no ssh"
            elif not t2_gpu(tg, offer["num_gpus"]):
                why = "B2 GPU load or throughput"
            elif not t3_upload(tg):
                why = "B3 upload"
            elif not t4_setup(tg, offer["num_gpus"]):
                why = "B4 set-up or speed"
        except Exception as e:  # an error of this script: the box is not usable without a check
            why = "error in the test: %r" % e
        if why:
            if not fail_box(st, iid, d0, tg, why, offer, t_rent):
                log("STOP: instance %s is not gone after the destroy" % iid)
                print("RESULT DESTROY_CHECK %s" % iid)
                return 6
            continue
        if not t5_tests(tg):
            fail_box(st, iid, d0, tg, "B5 E6 tests (code fault)", offer, t_rent)
            log("STOP: B5 failed, a code fault: report with a patch proposal")
            print("RESULT B5_FAIL")
            return 4
        ok = start_chain(tg, gate, mode)
        st["boxes"].append({"id": iid, "machine": offer.get("machine_id"), "host": offer.get("host_id"),
                            "why": "PASS", "ssh": list(tg), "chain_started": ok, "dph": offer["dph_total"],
                            "t_rent": t_rent, "q16_mode": mode})
        save_state(st)
        log("PASS box %s (machine %s), ssh %s:%d, chain %s"
            % (iid, offer.get("machine_id"), tg[0], tg[1], "started" if ok else "NOT started"))
        print("RESULT PASS %s %s %d %s" % (iid, tg[0], tg[1], ok))
        return 0 if ok else 3
    log("STOP: failed boxes cost %.3f $ (budget %.2f $)" % (st["spent"], BUDGET))
    print("RESULT BUDGET")
    return 2


if __name__ == "__main__":
    sys.exit(main())
