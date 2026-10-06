#!/usr/bin/env python
"""
Which large files on the instance can be deleted? (PLAN 4.2, decided 2026-10-04: delete large files
when they are no longer needed.) Dry run by default: prints the list and the sizes; --apply deletes.
Run it on the instance in the study folder. Deletion only after Chris's go for the printed list.

A file is deletable only if ALL its conditions hold:
  m_*.vti (snapshots)  the run is DONE, domains.json and dissmap.json exist and are newer than every
                       snapshot of the run, and both JSON files are listed in --local_md5 (fetched and
                       verified on the laptop)
  geometry.vti         as m_*.vti
  checkpoint.pt,       the run is DONE and its summary.json is listed in --local_md5 (resume is no
  ckpt_*.pt            longer needed); for a run with snapshots (snap_phases > 0 in config.json, e.g. all
                       Sobol runs) also: POSTPROC_DONE exists, domains.json holds the end state ("final")
                       and domains.json is listed in --local_md5 (audit 3, item 3)
  init.pt              NEVER by default (it is small, 40 MB, and other runs start from it: the alpha 0.01
                       runs of the Sobol subset and the 50 / 100 mT jobs of mesh300 use --init_from; a
                       job file that is written later cannot be checked now). Only with --include_init,
                       and then only if the run is DONE and every job in jobs/*.txt that uses it is DONE.
Everything else is never touched.

    python cleanup_instance.py --local_md5 local.md5            # dry run
    python cleanup_instance.py --local_md5 local.md5 --apply    # after the go
"""
import argparse
import json
import pathlib
import shlex

HERE = pathlib.Path(__file__).resolve().parent


def dependents(jobs_dir):
    """init.pt path (relative to runs/) -> names of the jobs that start from it."""
    dep = {}
    for jf in sorted(pathlib.Path(jobs_dir).glob("*.txt")):
        for line in jf.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            name, _, args = line.strip().partition(" ")
            a = shlex.split(args)
            if "--init_from" in a:
                src = a[a.index("--init_from") + 1]
                dep.setdefault(str(pathlib.PurePosixPath(src)), set()).add(name)
    return dep


def checkpoint_ok(d, local):
    """May checkpoint.pt / ckpt_*.pt of run folder d be deleted? (conditions in the module doc)"""
    if not (d / "DONE").exists() or ("%s/summary.json" % d.name) not in local:
        return False
    try:
        snap = json.loads((d / "config.json").read_text()).get("snap_phases", 0) > 0
    except (OSError, ValueError):
        return False
    if not snap:
        return True
    try:
        final = "final" in json.loads((d / "domains.json").read_text())
    except (OSError, ValueError):
        return False
    return (d / "POSTPROC_DONE").exists() and final and ("%s/domains.json" % d.name) in local


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=str(HERE / "runs"))
    ap.add_argument("--jobs", default=str(HERE / "jobs"))
    ap.add_argument("--local_md5", required=True, help="md5 list of the files on the laptop (paths relative to runs/)")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--include_init", action="store_true",
                    help="also delete init.pt (only if no pending job in jobs/*.txt uses it)")
    a = ap.parse_args()
    runs = pathlib.Path(a.runs)
    local = {line.split(None, 1)[1].strip() for line in open(a.local_md5) if line.strip()}
    dep = dependents(a.jobs)
    todo, keep = [], []
    for d in sorted(p for p in runs.iterdir() if p.is_dir()):
        done = (d / "DONE").exists()
        vtis = sorted(d.glob("*.vti"))
        if vtis:
            js = [d / "domains.json", d / "dissmap.json"]
            ok = done and all(j.exists() for j in js) and \
                all(j.stat().st_mtime > max(v.stat().st_mtime for v in vtis) for j in js) and \
                all("%s/%s" % (d.name, j.name) in local for j in js)
            (todo if ok else keep).extend((v, "snapshot", ok) for v in vtis)
        for ck in list(d.glob("checkpoint.pt")) + list(d.glob("ckpt_*.pt")):
            ok = checkpoint_ok(d, local)
            (todo if ok else keep).append((ck, "checkpoint", ok))
        ip = d / "init.pt"
        if ip.exists():
            users = dep.get("%s/init.pt" % d.name, set())
            pending = sorted(u for u in users if not (runs / u / "DONE").exists())
            ok = a.include_init and done and not pending
            why = "init (needed by %s)" % ", ".join(pending) if pending else "init (kept by default)"
            (todo if ok else keep).append((ip, why if not ok else "init", ok))
    size = lambda lst: sum(p.stat().st_size for p, _, _ in lst) / 1e9
    print("deletable: %d files, %.2f GB" % (len(todo), size(todo)))
    groups = {}
    for p, kind, _ in todo:
        groups.setdefault((p.parent.name, kind.split(" ")[0]), []).append(p)
    for (run, kind), ps in sorted(groups.items()):
        print("   %-24s %-10s %4d files  %.2f GB" % (run, kind, len(ps), sum(p.stat().st_size for p in ps) / 1e9))
    print("kept (conditions not met): %d files, %.2f GB" % (len(keep), size(keep)))
    for p, kind, _ in keep:
        if kind.startswith("init (needed"):
            print("   %s  %s" % (p.relative_to(runs), kind))
    n_init = sum(1 for _, k, _ in keep if k.startswith("init"))
    print("   of these init.pt: %d (never deleted without --include_init)" % n_init)
    if a.apply:
        for p, _, _ in todo:
            p.unlink()
        print("deleted %d files" % len(todo))
    else:
        print("dry run: nothing deleted (use --apply after the go)")


if __name__ == "__main__":
    main()
