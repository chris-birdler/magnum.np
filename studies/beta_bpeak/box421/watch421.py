#!/usr/bin/env python3
"""PLAN 4.21, local watcher of the box (runs on the laptop; the box chain runs without it). Every POLL seconds: the
state of the box (chain421.log, dispatcher status, the Bloch samples of the runs, the passes of the preparation, the
GPU and host memory, the end marker) and the cost. The script exits with one line "EVENT <kind> <text>" at the first
new event, so that the session gets a notification; start it again after the event is handled (the reported events
are kept in watch421_state.json).
Events (PLAN 4.21, "Early warnings and stops"):
  cost75     the cost of this plan reached 75 % of the limit (with a projection);
  costlimit  the cost reached the limit: stop, download, md5, destroy;
  credit     the credit is below the rest of the cost limit plus 1 $;
  bloch:<run> a Bloch sample in a cycle of a run of K1, K3 or K4 (decision D5);
  pass3:<run> a preparation needs its third pass (its second pass ended with a Bloch point);
  r4fail     the first-run check (R4) of a run failed (the chain starts no new jobs);
  failed:<n> a job failed (dispatcher status);
  gpumem     a GPU above 28 GiB (memory.used);
  hostmem    less than 8 GB of host memory available;
  dispatcher the dispatcher has ended or has disabled a GPU;
  dead       chain421.sh or dispatcher.sh does not run, and the chain has no end marker;
  end        the chain has ended (runs/CHAIN421_END);
  nossh      no ssh for 60 min;
  status     a status every STATUS_H hours.
Usage: watch421.py <instance id> <ssh host> <ssh port> <cost limit $> [<spent on failed boxes, $>]."""
import json, pathlib, re, subprocess, sys, time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import auto421 as A  # noqa: E402

STATE = HERE / "watch421_state.json"
POLL = 1200
STATUS_H = 3.0
B = "/root/magnum.np/studies/beta_bpeak/"
STATUS = "/root/magnum.np-disp/studies/beta_bpeak/status_dispatcher.json"
GPU_MAX_MIB = 28 * 1024
HOST_MIN_GB = 8


def credit():
    r = subprocess.run([A.W.V, "show", "user", "--raw"], capture_output=True, text=True, timeout=90)
    try:
        return float(json.loads(r.stdout)["credit"])
    except Exception:
        return None


def box_state(tg):
    cmd = ("tail -n 30 %sruns/chain421.log 2>/dev/null; echo '@@STATUS'; cat %s 2>/dev/null; echo '@@BLOCH'; "
           "for f in %sruns/K*_v16/bloch_[0-9][0-9]_amp[0-9][0-9].csv; do [ -f $f ] && "
           "echo \"$f $(awk -F, '!/^#/ && ($3 + $4) > 0' $f | wc -l)\"; done; echo '@@PASS'; "
           "for f in %sruns/K*_v16/stdout.log; do grep -q 'pass 1 (h0.*Bloch points +[0-9]* -[0-9]*' $f 2>/dev/null && "
           "grep 'pass 1 (h0' $f | grep -qv 'Bloch points +0 -0' && echo $f; done; echo '@@MEM'; "
           "nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits; "
           "if [ -r /sys/fs/cgroup/memory.max ]; then echo cg $(cat /sys/fs/cgroup/memory.max) $(cat /sys/fs/cgroup/memory.current); "
           "else free -b | awk '/Mem:/ {print \"cg\", $2, $3}'; fi; "
           "echo alive $(pgrep -fc '[b]ash chain421.sh') $(pgrep -fc '[b]ash /root/dispatcher.sh'); "
           "echo '@@END'; cat %sruns/CHAIN421_END 2>/dev/null; echo '@@DONE'; ls %sruns/*/DONE 2>/dev/null | wc -l"
           % (B, STATUS, B, B, B, B))
    rc, out = A.ssh(tg, cmd, timeout=120)
    if rc != 0 or "@@DONE" not in out:
        return None
    return parse_state(out)


def parse_state(out):
    """The box state from the output of the command of box_state."""
    log, rest = out.split("@@STATUS", 1)
    st, rest = rest.split("@@BLOCH", 1)
    bl, rest = rest.split("@@PASS", 1)
    ps, rest = rest.split("@@MEM", 1)
    mem, rest = rest.split("@@END", 1)
    end, ndone = rest.split("@@DONE", 1)
    try:
        status = json.loads(st)
    except Exception:
        status = {}
    bloch = {}
    for line in bl.strip().splitlines():
        m = re.match(r"\S*/(K\d_\d_v16)/(bloch_\S+\.csv) (\d+)", line)
        if m and int(m.group(3)) > 0:
            bloch.setdefault(m.group(1), {})[m.group(2)] = int(m.group(3))
    pass3 = sorted({re.sub(r".*/(K\d_\d_v16)/stdout.log", r"\1", x) for x in ps.split() if x.endswith("stdout.log")})
    gpu = {int(a): int(b) for a, b in re.findall(r"^(\d+),\s*(\d+)\s*$", mem, re.M)}
    cg = re.search(r"^cg (\S+) (\d+)", mem, re.M)                 # cgroup limit and use (bytes), else the host
    avail = None
    if cg and cg.group(1).isdigit():
        avail = (int(cg.group(1)) - int(cg.group(2))) / 1e9
    al = re.search(r"^alive (\d+) (\d+)", mem, re.M)
    return {"log": log.strip(), "status": status, "bloch": bloch, "pass3": pass3, "gpu_mib": gpu,
            "host_avail_gb": avail, "alive": (int(al.group(1)), int(al.group(2))) if al else None,
            "end": end.strip(), "n_done": int(ndone.strip() or 0)}


def main():
    iid, host, port, limit = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
    spent_failed = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    tg = (host, port)
    S = json.loads(STATE.read_text()) if STATE.exists() else {"reported": [], "last_status": 0.0, "last_ssh": time.time()}

    def event(kind, text):
        S["reported"].append(kind)
        STATE.write_text(json.dumps(S, indent=1))
        print("EVENT %s %s" % (kind, text), flush=True)
        sys.exit(0)

    while True:
        d = A.inst(iid) or {}
        if d.get("dph_total") and d.get("start_date"):          # keep the last known values if the call fails
            S["dph"], S["t0"] = float(d["dph_total"]), float(d["start_date"])
        dph = float(S.get("dph") or 0.0)
        t0 = float(S.get("t0") or time.time())
        hours = (time.time() - t0) / 3600.0
        cost = spent_failed + hours * dph
        cr = credit()
        bs = box_state(tg)
        now = time.time()
        if bs is None:
            if now - S["last_ssh"] > 3600 and "nossh" not in S["reported"]:
                event("nossh", "no ssh for %.0f min (instance %s, status %s)" % (
                    (now - S["last_ssh"]) / 60, iid, d.get("actual_status")))
        else:
            S["last_ssh"] = now
        info = "cost %.2f $ of %.0f $ (%.1f h at %.3f $/h, failed boxes %.2f $), credit %s" % (
            cost, limit, hours, dph, spent_failed, "%.2f $" % cr if cr is not None else "?")
        if cost >= limit and "costlimit" not in S["reported"]:
            event("costlimit", info)
        if cost >= 0.75 * limit and "cost75" not in S["reported"]:
            event("cost75", info + "; projection: see the run times in the chain log")
        if bs is not None:
            st = bs["status"]
            for n in sorted(bs["bloch"]):
                if n[:2] in ("K1", "K3", "K4") and "bloch:" + n not in S["reported"]:
                    event("bloch:" + n, "%s: Bloch samples %s (D5); %s" % (n, bs["bloch"][n], info))
            for n in bs["pass3"]:
                if "pass3:" + n not in S["reported"]:
                    event("pass3:" + n, "%s: the preparation needs its third pass; %s" % (n, info))
            if "R4 FAIL" in bs["log"] and "r4fail" not in S["reported"]:
                event("r4fail", "%s; %s" % ([x for x in bs["log"].splitlines() if "R4" in x][-1:], info))
            for n in sorted((st.get("failed") or {})):
                if "failed:" + n not in S["reported"]:
                    event("failed:" + n, "%s: %s; %s" % (n, st["failed"][n], info))
            if bs["gpu_mib"] and max(bs["gpu_mib"].values()) > GPU_MAX_MIB and "gpumem" not in S["reported"]:
                event("gpumem", "GPU memory used (MiB) %s; %s" % (bs["gpu_mib"], info))
            if bs["host_avail_gb"] is not None and bs["host_avail_gb"] < HOST_MIN_GB and "hostmem" not in S["reported"]:
                event("hostmem", "host memory available %.1f GB (cgroup); %s" % (bs["host_avail_gb"], info))
            if bs["alive"] and not bs["end"] and min(bs["alive"]) == 0 and "dead" not in S["reported"]:
                event("dead", "chain421.sh / dispatcher.sh processes %s, no end marker; %s" % (bs["alive"], info))
            if (st.get("final") or st.get("disabled_gpus")) and "dispatcher" not in S["reported"]:
                event("dispatcher", "final %s, disabled GPUs %s; %s" % (st.get("final"), st.get("disabled_gpus"), info))
            if bs["end"] and "end" not in S["reported"]:
                event("end", "%s; %s" % (bs["end"].replace("\n", " | "), info))
        if cr is not None and "credit" not in S["reported"] and cr < (limit - cost) + 1.0:
            event("credit", "credit %.2f $ below the rest of the cost limit plus 1 $; %s" % (cr, info))
        if now - S["last_status"] > STATUS_H * 3600:
            S["last_status"] = now
            STATE.write_text(json.dumps(S, indent=1))
            st = (bs or {}).get("status", {})
            print("EVENT status %s; running %s; done %s; waiting %s; GPU MiB %s; host avail %s GB; bloch %s; chain: %s"
                  % (info, sorted((st.get("running") or {})), (bs or {}).get("n_done"), st.get("waiting"),
                     (bs or {}).get("gpu_mib"), (bs or {}).get("host_avail_gb"), (bs or {}).get("bloch"),
                     ((bs or {}).get("log") or "").splitlines()[-1:]), flush=True)
            sys.exit(0)
        time.sleep(POLL)


if __name__ == "__main__":
    main()
