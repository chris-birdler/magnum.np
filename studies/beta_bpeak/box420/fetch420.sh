#!/bin/bash
# PLAN 4.20: download the files of the runs of a stage from the box and check them with md5 (from fetch416.sh).
# Stage B: the small files, the relaxed states (init.pt) and the field files (diag_fields.npz) of C8, D8, D8V1, D8V2,
# F8 (about 1 GB). Stage C: the small files and the relaxed states of the other runs (R0, HFP_s922, Dv, q 16, spare;
# about 0.8 GB), and the chain files. No checkpoints, no uploaded start states (D15, HFP_s921/923/924/925).
# Usage: fetch420.sh <ssh host> <ssh port> <B|C>. Local target: runs/plan420_box of the diag worktree; the md5 list of
# the box is kept there as MANIFEST_plan420_<stage>.md5. Exit code 0 only if each file agrees.
H=root@$1; P=$2; ST=$3
B=/root/magnum.np/studies/beta_bpeak/runs
L=/home/christoph/Documents/code/magnum.np-diag/studies/beta_bpeak/runs/plan420_box
T=$(mktemp -d)
[ -n "$1" ] && [ -n "$P" ] && { [ "$ST" = B ] || [ "$ST" = C ]; } || { echo "usage: fetch420.sh <host> <port> <B|C>"; exit 2; }
mkdir -p $L
if [ "$ST" = B ]; then SEL='C8_s* D8_s* D8V1_s* D8V2_s* F8_s*'; else
    SEL='R0_s* HFP_s922 DV1_s* DV2_s* D16_s* DV116_s* DV216_s* D_s925_*'; fi
SSHO="-o ConnectTimeout=20 -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null -o BatchMode=yes"
ssh $SSHO -p $P $H "cd $B && for r in $SEL; do [ -d \$r ] || continue; find \$r -maxdepth 1 -type f ! -name 'checkpoint*.pt' ! -name 'ckpt_*.pt' ! -name '*.vti' ! -name '*.tmp*' -print0 | xargs -0 -r md5sum; done; ls $( [ $ST = C ] && echo chain420.log ) SELECTED CHAIN420_END R0_FAIL plan420_stage*.json plan420_stage*.md 2>/dev/null | xargs -r md5sum" > $T/remote.md5
rc=$?; [ $rc -eq 0 ] || { echo "remote md5 failed rc $rc"; exit 1; }
cut -c35- $T/remote.md5 > $T/files.txt
rsync -a -e "ssh $SSHO -p $P" --files-from=$T/files.txt $H:$B/ $L/
rc=$?; [ $rc -eq 0 ] || { echo "rsync failed rc $rc"; exit 1; }
cd $L && md5sum -c --quiet $T/remote.md5; rc=$?
cp $T/remote.md5 $L/MANIFEST_plan420_$ST.md5
echo "stage $ST: files $(wc -l < $T/files.txt), $(du -sh $L | cut -f1) in $L, runs $(cut -d/ -f1 $T/files.txt | sort -u | wc -l), md5 check rc $rc"
exit $rc
