# HANDOVER: staircase speed patches, PLAN 14 (2026-10-24, about 14:20 UTC)

This note is for the Claude Code session that continues this work on another computer (agent handover, requested by
Chris). The record of the study stays PROTOKOLL.md (section 3.40). PLAN.md sections 14 and 14a hold the approved plan.
This note gives the state and the next step. It does not replace the PROTOKOLL.

## 1. Rules of Chris

- `~/.claude/CLAUDE.md` and the memory folder of the project hold the binding rules. Read both first.
- Chat in German, short. Documents, code and commits in English (ASD-STE100).
- Report a finding and its cause before a fix, then wait for the go. A code fix only with the word "patch" and a go.
- Commit and push after each finished step. No AI attribution in commits (memory no-ai-attribution).
- Destroy, release, merge and delete only with an explicit go. Before a destroy: download, md5, report, wait.
- Rent no box without a go. A go for "a new box" covers one box.
- Give each effect with its count n and its error type. Claim no effect from fewer than 3.
- PROTOKOLL.md is shared with other sessions: replace only up to the next heading. 3.40 is the last section.
- The GitHub repository chris-birdler/magnum.np is public. The paper repository and magnum.pi stay local.

## 2. Mandate

Chris (PLAN 14): all speed patches are tested, then the final speed gain is measured. His words: be safe and
conservative, the speed changes must not change the physics; choose the criteria; at the end, say which tests were
used and why they are sufficient and necessary.

## 3. State

- Done and PASS: T1, T2, T5, T6; T4 (K1) pairs A2, A1, A3 and AB off. Details: PROTOKOLL 3.40. The T1 and T2 results
  of PLAN 10, B5, B4, B1, A4 and AB on are for the CSR code; with the patch they must run again (section 4, step 9).
- Open: the B side of K1 (pairs PLAN 10, B5, B4, B1, A4, AB on), BA4 handling (window 48 has no determined point),
  then the speed measurement (PLAN 14, "Speed measurement").
- Finding (cause found): transient faults of the sparse product of the correction. The product is the CSR SpMV of
  `self._A @ m.reshape(-1)` in `DemagStaircaseCorrection.h` (magnumnp/field_terms/demag_staircase.py, from PLAN 10,
  01ed90c). torch uses cuSPARSE with CUSPARSE_SPMV_ALG_DEFAULT (from the torch source; the dispatch on the box is not
  verified).
  - Outlier test 3 (each try done twice in the real dynamics, 2 processes on the GPU): d538bde 5 events (H_d53 1, first
    tool, try count not saved; H2_d53 4 in 68 266 try evaluations); c63e862 (old pair sum with index_add_, no
    cuSPARSE) 0 in 37 920. Pooled with the pilot and the replay: d538bde 8 in 116 690, c63e862 0 in 61 760.
  - These 5 events are all in the correction; the 4 events of the refined tool are in the sparse product, each in 1 or
    2 entries. The 3 pilot events and the 1 replay event were not located. 3 of the 5 event rows have the length 513
    (0.47 % of the rows; probability about 1e-6 … 1e-5; after the bound for a pattern found after the fact, at most
    0.014 with row weights and 0.004 with entry weights).
  - 0 of 5 faults were accepted by RKF45, but a fault below about 1e-5 in a late stage would be accepted.
  - Stress test 3 (the product alone, back to back): 0 events in 3 408 200 products. Expected at least 17. Thus the
    fault needs the conditions of the real dynamics. A test of a fix must use the real dynamics.
  - All observed faults occurred with 2 or 3 processes on the GPU. ECC counters stay 0.
- Decision of Chris (2026-10-24, between 13:21 and 14:00 UTC, "go"): option 1, a patch that replaces the cuSPARSE
  product with an own sum per row in a fixed order, without cuSPARSE and without atomic additions. Verify it with
  outlier test 3 on the patched tree, then K1. Box 11 stays rented until the patch is tested. (Option 2 was: take PLAN
  10 out of the final state; option 3: accept with a bound, not recommended.) Recorded in PROTOKOLL 3.40.
- Nothing of the patch is built. I had read the code and had a draft design (section 4). The design is not reviewed.

## 4. Next step: the patch (draft design, not reviewed)

1. Branch `perf/staircase-det` from df0f811 (perf/staircase-speed). An unproven idea goes on its own branch.
2. Data layout. `pair_matrix` builds 3 × 3 symmetric blocks: block row t (target cell) has K_t blocks; each block has
   6 coefficients (SYM order) and 1 source cell. Start from the merged pairs: `_pair_matrix_cpu` adds the coefficients
   of duplicate (t, s) pairs into one block (sort, unique_consecutive, index_add_). Row length 3 K_t; K_t is 1 … 240,
   median 51 (geometry B, N 108, q 4). Draft: cut each block row into chunks of a fixed length C (for example 16).
   Pad the last chunk with a source index that points to an added zero row of m. Stage 1: gather m, `sym_matvec`, sum
   over the C blocks: partial sums [n_chunks, 3]. Stage 2: a padded table of the chunk indices of each block row;
   gather and sum: h of the rows. `h[rows] = ...` writes unique rows (no accumulation).
3. All operations are torch gather, element-wise operations and `sum(dim)`. These are deterministic and use no atomic
   additions. The order of the sum changes, thus the result changes at the round-off level only.
4. Estimates (not measured): the stored data about 1.4 GB against 2.6 GB of the CSR matrix (float64, geometry B,
   N 108). Eager mode adds temporary arrays of about 0.6 GB each (gathered m, products); process the chunks in
   blocks if needed. Speed references: the CSR product in float64 on geometry B, N 108 takes about 3.7 ms with the GPU
   alone and about 8 ms with 2 processes (stress test 3); in single precision on the main geometry at N 104: CSR 2.2 ms,
   old pair sum 20.8 ms (C4 table, PROTOKOLL 3.38). In df0f811 only llg.py uses `@torch.compile`; module functions
   with `@torch.compile` are in demag.py and rkf45.py of perf/timestep (c02f40b) and AB (7f755a0). A compiled module
   function can fuse the stages. Measure eager and compiled on box 11.
5. Before the build: let the advisor or a context-free agent review the design.
6. Local tests (CPU, float64): the new sum against the CSR product and against the old pair sum `NearCorrection.h_pair`
   (round-off level); bitwise equal on repetition; the existing staircase tests.
7. The tools hook the product by the attribute `_A`: `step_dump.wrap_terms` (`getattr(corr, "_A")`, `SparsePart`),
   `glitch_stress.py` and `csr_hist.py` (`corr._A.crow_indices()`). If the patch removes `_A`, the part "spmv"
   disappears without an error. Give the new sum its own hook (for example a method that the tools wrap), and make
   the tools fail loudly when they find neither.
8. Box 11: speed and bitwise repetition on the GPU against the CSR product.
9. Patched trees and repeated tests (scope for the go of Chris, see step 12). The file differs between the trees
   (d538bde, 01ed90c, c56c0a5, df0f811, BA4 b750b1d, AB 7f755a0): port the patch to each tree as its own commit on its
   own branch (for example perf/staircase-det-d53 …), and give each box tree a name with the suffix P (B_d53P …).
   T1 and T2 must run again for the pairs with the patch (PLAN 10, B5, B4, B1, A4, AB on); the speed measurement
   uses the patched AB.
10. Outlier test 3 on d538bde with the patch (B_d53P), in parallel with d538bde without the patch as the positive
    control (2 processes on the GPU, the condition of all observed faults). Job files: studies/staircase_correction/
    box11/h_redo.txt and h2_redo.txt (d538bde used `--k1_reps` 20 and 18, about 5 h). Draft pass rule (needs the go
    of Chris): the deterministic sum allows a bitwise check, so each try done twice must give bitwise equal results
    in the patched arm (the other terms were bitwise equal in glitch_test); verdict PASS when the patched arm has 0
    events in at least as many try evaluations as H2_d53 (68 266) and the control has at least 3 events in the same
    time; when the control has fewer than 3 events, extend both arms.
11. K1 with the patch: repeat the B side with the patch on each tree that had the CSR product: reference, base,
    control, 3 null runs, BA4 base, BA4 control, 3 BA4 null runs, 01ed90c, c56c0a5, df0f811, AB on (15 arms of
    10 … 15 min, 2 at a time, about 1.5 h). The arm B_c63 (c63e862, no sparse product) is done and stays valid.
12. Cost (estimate): outlier test about 5 h, K1 about 1.5 h, T1 and T2 about 1 h, speed tests about 1 h: about 9 h
    of box 11, about 1.5 $. Show steps 9 … 11 with their pass rules to Chris and get his go before the box runs
    (memory feedback-show-plan-before-run).
13. Then the K1 pairs, BA4 handling, the speed measurement and the final report (tests used, why sufficient and
    necessary, speed gain ± SE).

## 5. Box 11

- vast.ai instance 57396637, Tesla V100-SXM2 16 GB, 0.162 $/h. Access: see memory staircase-instance, or
  `vastai show instance 57396637 --raw`. Python: /opt/conda/bin/python3.11. Rented by this work; destroy only with go.
- No load since 2026-10-24 13:21 UTC (Chris approved to keep it until the patch is tested).
- Trees: /root/{B_c63,B_01e,B_c56,B_d53,B_df0,BA4,AB,A_31a,A_a1f,A_3d9,...}. Scripts: /root/s/ (k1.sh JOBFILE PAR,
  k1_eval.sh pilot|pairs, step_dump.py, glitch_stress.py, stress3.sh, find_try.py, m_end.py). Log of all runs:
  /root/progress.txt.
- K1 data: /root/k1/A/ and /root/k1/B/ (refB.pt 2.1 GB, arms about 1.6 GB each; not downloaded). Before a destroy:
  state the size, download what is needed, md5.
- Start long runs with `setsid nohup ... > log 2>&1 < /dev/null & disown`. Find a run's PID via /proc/PID/fd/1 (the
  process title has no script arguments). Kill by PID only, never `pkill -f`.

## 6. Local data and branches

- Raw data (gitignored, with md5 manifests): remote_20261023_box11/ (H/, csr/, stress3/, K1B/) and the older
  remote_* folders, 1.6 GB in total.
- Branches (remote "github"): feature/staircase-correction (PROTOKOLL, PLAN, study scripts), perf/staircase-speed
  df0f811, perf/staircase-ba4 b750b1d, perf/speed-all 7f755a0 (AB), perf/timestep c02f40b, perf/a1-fixed a47492b.
- Other open points (not urgent): 2-thread CPU non-reproducibility (not examined); option 2 parked (memory
  staircase-option2-parked); paper draft (memory paper-pbc-staircase).
