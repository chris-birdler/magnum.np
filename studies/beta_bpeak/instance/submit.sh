#!/bin/bash
# Submit a job file to the box dispatcher (dispatcher.sh): submit.sh <priority> <tag> <job file>. A smaller priority
# runs first. The dispatcher writes spool/<priority>_<tag>.jobs.done when all jobs of the file are finished.
B=/root/magnum.np/studies/beta_bpeak
[ $# -eq 3 ] || { echo "usage: submit.sh <priority> <tag> <job file>"; exit 1; }
f=$B/spool/$1_$2.jobs
[ -e $f ] && { echo "$f exists: use another tag"; exit 1; }
mkdir -p $B/spool && cp "$3" $f.tmp && mv $f.tmp $f && echo "submitted $f ($(grep -c . $f) jobs)"
