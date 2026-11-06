#!/usr/bin/env python3
"""PLAN 4.21 (from watch_rent420.py of PLAN 4.20): V100 32 GB, exactly 4 GPUs, at most 0.10 $ per GPU-h, CUDA 12.6 or
newer, not deverified, host memory at least 128 GB, disk 160 GB (PLAN 4.21, "Box"); not the machines 147583, 152226,
152674, 154925 and not the host 381795.
Watch vast.ai V100 offers; rent the first one that meets the rule of Chris (2026-10-27, raised the same day for
PLAN 4.16): 4x V100 at 0.08 $ per GPU-h or less (0.321 $/h, the price of machine 147583); 16 GB allowed
(Chris 2026-10-28, option 2: the speed reference box 11 was a V100-SXM2 16 GB; 156^3 with q 8 needed 10.4 GiB).
rent_one() rents at most one instance (auto421.py adds the box tests). Log: watch_rent421.log next to this file."""
import json, pathlib, subprocess, sys, time

V = "/home/christoph/Documents/.venv/bin/vastai"
LOG = pathlib.Path(__file__).with_suffix(".log")  # watch_rent421.log
QUERY = ("gpu_name in [Tesla_V100,V100_SXM2,Tesla_V100S] num_gpus=4 gpu_ram>=31 rentable=true verified=any "
         "reliability>=0.98 disk_space>=200 cuda_max_good>=12.6 inet_down>=100")
LIMIT = 0.1005        # $ per GPU-h (dph_total of the search, default storage); 0.10 (Chris 2026-11-03; before 0.094)
RAM_MIN_GB = 128.0    # host memory in GB (PLAN 4.21: one q-16 set-up of 34.5 GiB at a time beside 3 runs); cpu_ram in MiB
IMAGE = "pytorch/pytorch:2.7.1-cuda12.6-cudnn9-devel"
DISK = 160            # GB (PLAN 4.21: snapshots up to about 9 GB per run until the post-processing, 4 runs at a time)
LABEL = "beta-plan421"
MAX_HOURS = 72
POLL = 120
EXCLUDED = {147583, 152226, 152674, 154925}  # 147583: could not pull IMAGE (2026-10-27); 152226: foreign GPU load; 152674: GPU 0 thermal slowdown (2026-10-31); 154925: speed check failed (2026-11-01)
EXCLUDED_HOSTS = {381795}  # all machines of this host (one public IP, 166.113.48.43): same network (Chris, 2026-10-28)


def log(msg):
    with LOG.open("a") as f:
        f.write(time.strftime("%Y-%m-%d %H:%M:%S UTC ", time.gmtime()) + msg + "\n")


def search():
    r = subprocess.run([V, "search", "offers", QUERY, "-o", "dph_total", "--raw"], capture_output=True, text=True,
                       timeout=90)
    return json.loads(r.stdout) if r.returncode == 0 and r.stdout.strip().startswith("[") else None


def rent_one(excluded=(), max_hours=MAX_HOURS):
    """Poll until one offer is rented; excluded: more machine ids. Returns (instance id, offer) or None."""
    t0 = time.time()
    excl = EXCLUDED | set(excluded)
    log("start: limit %.3f $/GPU-h, excluded machines %s, hosts %s, query %s" % (LIMIT, sorted(excl), sorted(EXCLUDED_HOSTS), QUERY))
    n = 0
    while time.time() - t0 < max_hours * 3600:
        try:
            offers = search()
        except Exception as e:  # network or CLI error: try again later
            offers = None
            log("search error: %r" % e)
        if offers is not None:
            ok = [x for x in offers if x["dph_total"] / x["num_gpus"] < LIMIT
                  and x.get("cpu_ram", 0) * 2**20 / 1e9 >= RAM_MIN_GB
                  and x.get("verification") != "deverified"]   # a host that lost its verification: not rented
            ok = [x for x in ok if x.get("machine_id") not in excl and x.get("host_id") not in EXCLUDED_HOSTS]
            ok.sort(key=lambda x: (x["dph_total"] / x["num_gpus"], -x["num_gpus"]))
            best = min((x["dph_total"] / x["num_gpus"] for x in offers), default=float("nan"))
            n += 1
            if n % 15 == 1 or ok:
                log("poll %d: %d offers, best %.3f $/GPU-h, %d below the limit" % (n, len(offers), best, len(ok)))
            for x in ok:
                desc = "offer %d: %dx %s %.0f GB, %.3f $/h (%.3f $/GPU-h), RAM %.0f GiB, rel %.3f, %s, %s" % (
                    x["id"], x["num_gpus"], x["gpu_name"], x["gpu_ram"] / 1024, x["dph_total"],
                    x["dph_total"] / x["num_gpus"], x["cpu_ram"] / 1024, x.get("reliability2", 0),
                    x.get("geolocation", "?"), x.get("verification", "?"))
                log("try " + desc)
                r = subprocess.run([V, "create", "instance", str(x["id"]), "--image", IMAGE, "--disk", str(DISK),
                                    "--label", LABEL, "--ssh", "--direct", "--cancel-unavail", "--raw"],
                                   capture_output=True, text=True, timeout=120)
                log("create rc %d: %s %s" % (r.returncode, r.stdout.strip()[:400], r.stderr.strip()[:200]))
                try:
                    res = json.loads(r.stdout)
                except Exception:
                    res = {}
                if res.get("success") and res.get("new_contract"):
                    log("RENTED instance %s (%s)" % (res["new_contract"], desc))
                    return res["new_contract"], x
        time.sleep(POLL)
    log("stop: no offer below the limit in %d h" % max_hours)
    return None


def main():
    r = rent_one()
    if r is None:
        print("NO OFFER in %d h" % MAX_HOURS)
        return 1
    print("RENTED instance %s" % r[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
