# Decisions taken overnight 2026-10-19 / 20 (for the review of Chris)

Rule of Chris (2026-10-19, evening): if a problem comes up while he is away, write the options with a recommendation,
use the recommendation at once, and record each decision here and in the PROTOKOLL. Times UTC.

| time | problem | options | recommendation (used) | result |
|---|---|---|---|---|
| 22:01 | The dispatcher of the instance runs the code eba2dfd; the static-saturation jobs need the new code ebcea74. The K3_2 history runs still run. | (1) switch the dispatcher worktree to ebcea74 while it runs (as on 2026-10-19 10:14); the running jobs keep their loaded code, and the modules they import later did not change (checked by git diff); (2) stop the dispatcher with spool/EXIT, wait for K3_2 (about 30 min, 2 GPUs idle), restart. | (1) | Done; /root/disp_commit = ebcea74; the 8 jobs submitted at 22:01 (spool 25_satstatic); K4_1_nx_sad, _sbd, _sa200 started at once. **Confirm.** |
| 22:10 | The reading of the history test (7.123) was written one-sided for K4: "ln (w_l(ad) / w_l(nx)) > ln 1.5". At K3 the option C loss is lower than the nx loss. | (1) apply the rule literally: K3 "not determined" in all 3; (2) apply it two-sided, \|ln ratio\| > ln 1.5, because the purpose is "the two stable states differ by more than a factor 1.5". | (2) | K3_0 not determined (× 1.06), K3_1 and K3_2 two stable states (× 0.40, × 0.62); overall "two stable states" (2 of 3). Both readings are in PROTOKOLL 7.127. **Confirm.** |
