#!/usr/bin/env python3
"""PLAN 4.21 Part 1: the box tests B1 ... B5 and the chain start of auto420.py for a box that is already rented
(auto420.py was stopped after the rent). Difference to auto420.py (decision of Chris 2026-11-07, for this box only):
a fail of the GPU throughput (B2) or of the speed check (B4) does not destroy the box; the script writes HOLD with the
numbers and ends, and Chris decides. The other fails (B1 ssh, B2 GPU count or foreign load, B3 upload, B4 set-up, B5)
destroy the box as in auto420.py.
Usage: resume420.py <instance id> <commit> <gate A variants> <ram GB> [accept]
  without "accept": B1 ... B4, then (all pass) B5 and the chain start;
  with "accept": Chris accepted the B2 or B4 result: B5 and the chain start (the set-up of B4 is done).
Log: auto420.log (the functions of auto420.py)."""
import sys, time

import auto420 as A


def main():
    iid, A.COMMIT, gate, ram_gb = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    accept = len(sys.argv) > 5 and sys.argv[5] == "accept"
    mode = "setup" if ram_gb >= 128.0 else "serial"
    st = A.load_state()
    d0 = A.inst(iid)
    offer = {k: d0.get(k) for k in ("machine_id", "host_id", "dph_total", "num_gpus")}
    t_rent = float(d0.get("start_date") or time.time())
    n = int(offer["num_gpus"])
    A.log("resume420 %s: box %s (machine %s), %d GPUs, q16 mode %s, commit %s, gate %s"
          % ("accept" if accept else "tests", iid, offer["machine_id"], n, mode, A.COMMIT, gate))
    tg, why, hold = None, None, None
    try:
        tg = A.t1_ssh(iid)
        if tg is None:
            why = "B1 no ssh"
        elif not accept:
            mark = A.LOG.stat().st_size
            if not A.t2_gpu(tg, n):
                new = A.LOG.read_bytes()[mark:].decode()
                if "T2 FAIL: throughput" in new:
                    hold = "B2 throughput"
                else:
                    why = "B2 GPU load"
            elif not A.t3_upload(tg):
                why = "B3 upload"
            else:
                mark = A.LOG.stat().st_size
                if not A.t4_setup(tg, n):
                    new = A.LOG.read_bytes()[mark:].decode()
                    if new.count("T4 SPEED gpu") == n and "failed" not in new:
                        hold = "B4 speed"     # the set-up ended and each GPU has a speed value
                    else:
                        why = "B4 set-up"
    except Exception as e:
        why = "error in the test: %r" % e
    if hold:
        _, out = A.ssh(tg, "nvidia-smi --query-gpu=index,clocks.sm,clocks.max.sm,temperature.gpu,power.draw,"
                           "power.limit,clocks_event_reasons.active --format=csv,noheader", timeout=60)
        for line in out.strip().splitlines():
            A.log("HOLD nvidia-smi " + line)
        A.log("HOLD box %s (%s): the box stays rented without jobs; decision of Chris" % (iid, hold))
        print("RESULT HOLD %s %s %d" % (iid, tg[0], tg[1]))
        return 7
    if why:
        A.fail_box(st, iid, d0, tg, why, offer, t_rent)
        print("RESULT FAIL %s %s" % (iid, why))
        return 1
    if not A.t5_tests(tg):
        A.fail_box(st, iid, d0, tg, "B5 E6 tests (code fault)", offer, t_rent)
        A.log("STOP: B5 failed, a code fault: report with a patch proposal")
        print("RESULT B5_FAIL")
        return 4
    ok = A.start_chain(tg, gate, mode)
    st["boxes"].append({"id": iid, "machine": offer["machine_id"], "host": offer["host_id"], "why": "PASS",
                        "ssh": list(tg), "chain_started": ok, "dph": offer["dph_total"], "t_rent": t_rent,
                        "q16_mode": mode, "accepted_by_chris": accept})
    A.save_state(st)
    A.log("PASS box %s (machine %s), ssh %s:%d, chain %s" % (iid, offer["machine_id"], tg[0], tg[1],
                                                          "started" if ok else "NOT started"))
    print("RESULT PASS %s %s %d %s" % (iid, tg[0], tg[1], ok))
    return 0 if ok else 3


if __name__ == "__main__":
    sys.exit(main())
