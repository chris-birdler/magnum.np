"""Compare two benchmark result directories.

Usage: python compare.py results/base results/opt
"""
import json
import glob
import os
import sys

# metric -> True if higher is better
METRICS = {
    "h_ms": False,
    "init_s": False,
    "steps_per_s": True,
    "iters_per_s": True,
    "cuda_peak_mb": False,
    "ru_maxrss_mb": False,
}
REGRESSION_THRESHOLD = 3. # percent


def load(d):
    records = {}
    for path in sorted(glob.glob(os.path.join(d, "*.json"))):
        for record in json.load(open(path)):
            records[record["key"]] = record
    return records


def main():
    base_dir, opt_dir = sys.argv[1], sys.argv[2]
    base, opt = load(base_dir), load(opt_dir)

    regressions = []
    header = "%-36s %-14s %12s %12s %9s" % ("key", "metric", "base", "opt", "delta")
    print(header)
    print("-" * len(header))
    for key in sorted(base):
        if key not in opt:
            print("%-36s missing in %s" % (key, opt_dir))
            continue
        for metric, higher_is_better in METRICS.items():
            if metric not in base[key] or metric not in opt[key]:
                continue
            b, o = base[key][metric], opt[key][metric]
            if b == 0:
                continue
            delta = (o - b) / b * 100.
            improvement = delta if higher_is_better else -delta
            flag = ""
            if improvement < -REGRESSION_THRESHOLD:
                flag = "  << REGRESSION"
                regressions.append((key, metric, delta))
            print("%-36s %-14s %12g %12g %+8.1f%%%s" % (key, metric, b, o, delta, flag))

    print()
    if regressions:
        print("%d regression(s) worse than %g%%:" % (len(regressions), REGRESSION_THRESHOLD))
        for key, metric, delta in regressions:
            print("  %s %s (%+.1f%%)" % (key, metric, delta))
        sys.exit(1)
    print("no regressions worse than %g%%" % REGRESSION_THRESHOLD)


if __name__ == "__main__":
    main()
