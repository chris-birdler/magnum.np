#!/usr/bin/env python3
"""PLAN 4.21 Part 1: wait for an offer under the rent rules of watch_rent420.py (read only: search and credit, no rent).
Every 20 min one log line (box420/wait421p1.log): the offers, the eligible offers and the credit. The script ends with
one line "EVENT offer ..." when an offer meets the rules, or "EVENT timeout" after the given hours. At the offer the
session starts auto420.py; a low credit does not stop the rent (rule of Chris 2026-11-07), the session tells Chris the
credit. Usage: wait421p1.py [hours, default 24]"""
import json
import subprocess
import sys
import time

import watch_rent420 as W

POLL = 1200
LOG = W.LOG.with_name("wait421p1.log")


def credit():
    r = subprocess.run([W.V, "show", "user", "--raw"], capture_output=True, text=True, timeout=90)
    try:
        return float(json.loads(r.stdout)["credit"])
    except Exception:
        return None


def log(msg):
    with LOG.open("a") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S UTC ", time.gmtime()) + msg + "\n")


if __name__ == "__main__":
    hours = float(sys.argv[1]) if len(sys.argv) > 1 else 24.0
    t0 = time.time()
    log("start: %.0f h, query %s, limit %.4f $/GPU-h, RAM %.0f GB" % (hours, W.QUERY, W.LIMIT, W.RAM_MIN_GB))
    while time.time() - t0 < hours * 3600:
        try:
            offers = W.search()
        except Exception as e:
            offers = None
            log("search error: %r" % e)
        cr = credit()
        if offers is not None:
            ok = [x for x in offers if x["dph_total"] / x["num_gpus"] < W.LIMIT
                  and x.get("cpu_ram", 0) * 2**20 / 1e9 >= W.RAM_MIN_GB and x.get("verification") != "deverified"
                  and x.get("machine_id") not in W.EXCLUDED and x.get("host_id") not in W.EXCLUDED_HOSTS]
            ok.sort(key=lambda x: (x["dph_total"] / x["num_gpus"], -x["num_gpus"]))
            log("%d offers, %d eligible, credit %s" % (len(offers), len(ok), cr))
            if ok:
                x = ok[0]
                print("EVENT offer machine %s host %s %d GPUs %.4f $/GPU-h %.0f GB %s, credit %s" % (
                    x.get("machine_id"), x.get("host_id"), x["num_gpus"], x["dph_total"] / x["num_gpus"],
                    x.get("cpu_ram", 0) * 2**20 / 1e9, x.get("verification"), cr), flush=True)
                log("EVENT offer machine %s" % x.get("machine_id"))
                sys.exit(0)
        time.sleep(POLL)
    print("EVENT timeout %.0f h without an offer, credit %s" % (hours, credit()), flush=True)
    log("EVENT timeout")
