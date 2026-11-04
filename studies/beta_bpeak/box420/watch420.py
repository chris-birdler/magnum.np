#!/usr/bin/env python3
"""PLAN 4.20, local watcher of the box (runs on the laptop; the box chain runs without it). Every POLL seconds: the
state of the box (chain420.log, dispatcher status, the relaxation of the F8 jobs, the end markers) and the cost. The
script exits with one line "EVENT <kind> <text>" at the first new event, so that the session gets a notification;
start it again after the event is handled (the reported events are kept in watch420_state.json).
Events (PLAN 4.20, gates 4 ... 8):
  cost75     the cost of this plan reached 75 % of the limit (with a projection);
  costlimit  the cost reached the limit: stop, download, md5, destroy;
  credit     the credit is below the rest of the cost limit plus 1 $;
  f8warn     an F8 job at 20000 relaxation iterations without a stop (early warning);
  failed     a job failed (dispatcher status);
  selected   the selection of stage B is written (runs/SELECTED);
  end        the chain has ended (runs/CHAIN420_END);
  nossh      no ssh for 60 min;
  status     a status every STATUS_H hours.
Usage: watch420.py <instance id> <ssh host> <ssh port> [<spent on failed boxes, $>]."""
import json, pathlib, re, subprocess, sys, time

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import auto420 as A  # noqa: E402

STATE = HERE / "watch420_state.json"
LIMIT = 16.0                 # $ (decision of Chris, 2026-11-04)
POLL = 1200
STATUS_H = 3.0
B = "/root/magnum.np/studies/beta_bpeak/"
STATUS = "/root/magnum.np-disp/studies/beta_bpeak/status_dispatcher.json"
F8_WARN = 20000


def credit():
    r = subprocess.run([A.W.V, "show", "user", "--raw"], capture_output=True, text=True, timeout=90)
    try:
        return float(json.loads(r.stdout)["credit"])
    except Exception:
        return None


def box_state(tg):
    cmd = ("cat %sruns/chain420.log 2>/dev/null | tail -n 30; echo '@@STATUS'; cat %s 2>/dev/null; echo '@@F8'; "
           "for f in %sruns/F8_s94?/stdout.log; do [ -f $f ] && echo \"$f $(grep -c '^\\[relax\\] it' $f) "
           "$(grep -c 'virgin state: converged' $f)\"; done; echo '@@END'; "
           "cat %sruns/SELECTED %sruns/CHAIN420_END 2>/dev/null; echo '@@DONE'; ls %sruns/*/DONE 2>/dev/null | wc -l"
           % (B, STATUS, B, B, B, B))
    rc, out = A.ssh(tg, cmd, timeout=120)
    if rc != 0 or "@@DONE" not in out:
        return None
    return parse_state(out)


def parse_state(out):
    """The box state from the output of the command of box_state."""
    log, rest = out.split("@@STATUS", 1)
    st, rest = rest.split("@@F8", 1)
    f8, rest = rest.split("@@END", 1)
    end, ndone = rest.split("@@DONE", 1)
    try:
        status = json.loads(st)
    except Exception:
        status = {}
    f8it = {}
    for line in f8.strip().splitlines():
        m = re.match(r"\S*/(F8_s\d+)/stdout.log (\d+) (\d+)", line)
        if m:
            # chunks of 500 iterations (--relax_chunk 500); the stop line of the relaxation
            f8it[m.group(1)] = (int(m.group(2)) * 500, int(m.group(3)) > 0)
    return {"log": log.strip(), "status": status, "f8_iterations": f8it, "end": end.strip(),
            "n_done": int(ndone.strip() or 0)}


def main():
    iid, host, port = sys.argv[1], sys.argv[2], int(sys.argv[3])
    spent_failed = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    tg = (host, port)
    S = json.loads(STATE.read_text()) if STATE.exists() else {"reported": [], "last_status": 0.0, "last_ssh": time.time()}

    def event(kind, text):
        S["reported"].append(kind)
        STATE.write_text(json.dumps(S, indent=1))
        print("EVENT %s %s" % (kind, text), flush=True)
        sys.exit(0)

    while True:
        d = A.inst(iid) or {}
        dph = float(d.get("dph_total") or 0.0)
        t0 = float(d.get("start_date") or time.time())
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
            cost, LIMIT, hours, dph, spent_failed, "%.2f $" % cr if cr is not None else "?")
        if cost >= LIMIT and "costlimit" not in S["reported"]:
            event("costlimit", info)
        if cost >= 0.75 * LIMIT and "cost75" not in S["reported"]:
            event("cost75", info)
        if bs is not None:
            st = bs["status"]
            for n in sorted((st.get("failed") or {})):
                if "failed:" + n not in S["reported"]:
                    event("failed:" + n, "%s: %s; %s" % (n, st["failed"][n], info))
            for n, (it, stopped) in sorted(bs["f8_iterations"].items()):
                if it >= F8_WARN and not stopped and "f8warn:" + n not in S["reported"]:
                    event("f8warn:" + n, "%s at %d relaxation iterations without a stop; %s" % (n, it, info))
            if bs["end"]:
                first = bs["end"].splitlines()[0]
                if first in ("V1", "V2", "NONE") and "selected" not in S["reported"]:
                    event("selected", "selection %s; %s" % (first, info))
                if "CHAIN420_END" not in S["reported"] and any(x.startswith(("stage C done", "no stage C"))
                                                                 for x in bs["end"].splitlines()):
                    S["reported"].append("CHAIN420_END")
                    event("end", "%s; %s" % (bs["end"].replace("\n", " | "), info))
        if cr is not None and "credit" not in S["reported"] and cr < (LIMIT - cost) + 1.0:
            # the remaining jobs can cost up to the limit (gates 8)
            event("credit", "credit %.2f $ below the rest of the cost limit plus 1 $; %s" % (cr, info))
        if now - S["last_status"] > STATUS_H * 3600:
            S["last_status"] = now
            STATE.write_text(json.dumps(S, indent=1))
            st = (bs or {}).get("status", {})
            print("EVENT status %s; running %s; done %s; waiting %s; F8 iterations %s; chain: %s" % (
                info, sorted((st.get("running") or {})), (bs or {}).get("n_done"), st.get("waiting"),
                (bs or {}).get("f8_iterations"), ((bs or {}).get("log") or "").splitlines()[-1:]), flush=True)
            sys.exit(0)
        time.sleep(POLL)


if __name__ == "__main__":
    main()
