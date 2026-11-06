#!/bin/bash
# PLAN 4.21: download the files of the 12 runs from the box and check them with md5 (from fetch420.sh). Per run: the
# small files (config, summary, logs, samples, bloch, E3 files, geometry.npz, the results of the post-processing and of
# the checks), init.pt and the 7 stage states ckpt_*.pt (about 3 ... 4 GB in total); not checkpoint.pt (the end state
# is the last stage state), no snapshots (*.vti) and no accumulator files (acc_*.npz; the post-processing deletes them;
# a run whose post-processing failed keeps them on the box, and the report names it). Also the chain files.
# Usage: fetch421.sh <ssh host> <ssh port> [kept]. Local target: runs/plan421_box of the diag worktree; the md5 list of
# the box is kept there as MANIFEST_plan421.md5. Exit code 0 only if each file agrees.
# "kept": only the snapshots and the accumulator files of the finished runs without POSTPROC_DONE (an R4 FAIL or a
# failed evaluation, chain421.sh); up to about 9 GB per run. Run it first with "kept size" to print the size (PLAN 4.21:
# state the size before a download above 1 GB); md5 list MANIFEST_plan421_kept.md5.
H=root@$1; P=$2
if [ "$3" = kept ]; then
    B=/root/magnum.np/studies/beta_bpeak/runs
    L=/home/christoph/Documents/code/magnum.np-diag/studies/beta_bpeak/runs/plan421_box
    T=$(mktemp -d)
    SSHO="-o ConnectTimeout=20 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o BatchMode=yes"
    SEL="cd $B && for r in K*_v16; do [ -f \$r/DONE ] && [ ! -f \$r/POSTPROC_DONE ] && ls \$r/*.vti \$r/acc_*.npz 2>/dev/null; done"
    if [ "$4" = size ]; then
        ssh $SSHO -p $P $H "$SEL | xargs -r du -cb | tail -1"; exit $?
    fi
    ssh $SSHO -p $P $H "$SEL | xargs -r md5sum" > $T/remote.md5 || { echo "remote md5 failed"; exit 1; }
    cut -c35- $T/remote.md5 > $T/files.txt
    mkdir -p $L
    rsync -a -e "ssh $SSHO -p $P" --files-from=$T/files.txt $H:$B/ $L/ || { echo "rsync failed"; exit 1; }
    cd $L && md5sum -c --quiet $T/remote.md5; rc=$?
    cp $T/remote.md5 $L/MANIFEST_plan421_kept.md5
    echo "kept files $(wc -l < $T/files.txt), md5 check rc $rc"
    exit $rc
fi
B=/root/magnum.np/studies/beta_bpeak/runs
L=/home/christoph/Documents/code/magnum.np-diag/studies/beta_bpeak/runs/plan421_box
T=$(mktemp -d)
[ -n "$1" ] && [ -n "$P" ] || { echo "usage: fetch421.sh <host> <port>"; exit 2; }
mkdir -p $L
SSHO="-o ConnectTimeout=20 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o BatchMode=yes"
ssh $SSHO -p $P $H "cd $B && for r in K*_v16; do [ -d \$r ] || continue; find \$r -maxdepth 1 -type f ! -name 'checkpoint.pt' ! -name '*.vti' ! -name 'acc_*.npz' ! -name '*.tmp*' -print0 | xargs -0 -r md5sum; done; ls chain421.log postproc421.log CHAIN421_END STOP_NEW_JOBS 2>/dev/null | xargs -r md5sum" > $T/remote.md5
rc=$?; [ $rc -eq 0 ] || { echo "remote md5 failed rc $rc"; exit 1; }
cut -c35- $T/remote.md5 > $T/files.txt
rsync -a -e "ssh $SSHO -p $P" --files-from=$T/files.txt $H:$B/ $L/
rc=$?; [ $rc -eq 0 ] || { echo "rsync failed rc $rc"; exit 1; }
cd $L && md5sum -c --quiet $T/remote.md5; rc=$?
cp $T/remote.md5 $L/MANIFEST_plan421.md5
echo "files $(wc -l < $T/files.txt), $(du -sh $L | cut -f1) in $L, runs $(cut -d/ -f1 $T/files.txt | sort -u | wc -l), md5 check rc $rc"
exit $rc
