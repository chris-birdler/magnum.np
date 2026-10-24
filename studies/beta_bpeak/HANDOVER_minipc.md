# Handover: β(B_peak) corner study, laptop → mini PC (2026-10-24)

This note is for the next Claude Code session on the mini PC. It covers only the work of the laptop chat 6d24d835
(start folder `~/Documents/code/magnum.np/powder`): the β(B_peak) corner study (PLAN 4.15 phase 1), PROTOKOLL
7.138 … 7.141, the evaluation scripts on main, and the deck `beta_powder_sobol`. It does not cover the staircase
session (its own note: `studies/staircase_correction/HANDOVER.md`), the deck `via_3d_periodic` or the migration
itself (`~/Documents/umzug_minipc/`). The scientific record is `PROTOKOLL.md`. This note does not replace it and
contains no new results.

Terms in this note:
- "laptop chat": the chat 6d24d835 on the laptop that wrote this note.
- "instance": a rented vast.ai machine. The memory files also call it "box".
- "corner difference": β(Ki) − β(K1). PROTOKOLL 7.141 and §0 call it "lever against K1". The deck uses "lever" only
  for the Holm class of the main fit.
- "code part" and "preparation part": the share of a change that comes from the new code and from the new
  preparation (PROTOKOLL 7.141). The deck calls them "new code" and "new start".

## 1. Start and read first, in this order

Start Claude Code in `~/Documents/code/magnum.np/powder`. This folder loads the project memory below.

1. `~/.claude/CLAUDE.md` (working profile of Chris).
2. Project memory `~/.claude/projects/-home-christoph-Documents-code-magnum-np/memory/`: `MEMORY.md`, then
   `beta-bpeak-study.md`, `branch-layout.md`, `no-ai-attribution.md`, `feedback-short-answers.md`,
   `feedback-always-german.md`, `feedback-do-remote-ops-yourself.md`, `protokoll-shared-edits.md`.
3. PROTOKOLL (path in section 2): the status notes at the top of §0 and §12, then 7.138 … 7.141.
4. PLAN.md of the study, section 4.15 (corner study and new main study; D3 is open).

## 2. Where everything is

| Item | Path | Branch | HEAD at hand-over | Remote |
|---|---|---|---|---|
| Study record (PROTOKOLL, PLAN, results) | `~/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak/` (below: `<study>`) | `study/beta-bpeak-pinning` | the commit of this note (parent 1639f6a) | `github` (pushed) |
| Code (run and evaluation scripts, tests) | `~/Documents/code/magnum.np-main/studies/beta_bpeak/` | `main` | 4c6b7fe | `github` (pushed) |
| Deck | `~/Documents/slide_decks/beta_powder_sobol/` (`make_deck.py`, `PLAN.md`, `.pptx`, `.pdf`, `figures/`) | `main` | 90047d3 | `origin` (pushed) |
| Copy of the PROTOKOLL | `~/Documents/code/magnum.np/powder/PROTOKOLL.md` (not tracked) | – | identical to 1639f6a | – |
| Raw runs | `<study>/runs/` (26 GB, 915 entries, md5 lists `MANIFEST_*.md5` and `md5_*.txt`) | not tracked | – | local only |

Rules for these places:
- Edit and run code only in `magnum.np-main`. The study worktree has old copies of `run_loops.py` and
  `test_study.py` (code state 515b260) and no `eval_corners_compare.py`. Use the study worktree only for the
  PROTOKOLL, PLAN.md and `results/`.
- Do not commit in the checkout `~/Documents/code/magnum.np`. It is on `fix/preexisting-issues`, which Chris keeps
  out of main (`branch-layout.md`). Its folder `powder/` is not tracked.
- The PROTOKOLL copy in `powder/` is read-only. The laptop chat copied the study file to it after each PROTOKOLL
  commit. Continue this practice. Never edit the copy.
- In the magnum.np repository, `origin` is the gitlab upstream of the magnum.np team. Never push there. Push to
  `github`.
- `<study>/audit_20261012/sobol_audit_boot.npz` is untracked. It is generated output. Do not commit it.

Tools:
- `~/Documents/.venv/bin/python` runs the evaluation scripts, `test_study.py` and the deck build. This venv has no
  torch and no magnumnp. Thus `run_loops.py` and the magnum.np test suite do not run with it. To build a torch venv,
  see `magnumnp-testing.md`.
- vast.ai CLI: `~/Documents/.venv/bin/vastai` (use `show instances-v1`; `show instances` does not work).
- PDF export: `/snap/bin/libreoffice --headless --convert-to pdf` (`soffice` is not on the PATH).

## 3. State at hand-over

Done and pushed (2026-10-23):
- PLAN 4.15 phase 1 is complete (PROTOKOLL 7.138 … 7.140, study commits 08768aa and 212af0a).
- PROTOKOLL 7.141: comparison of the corner study with the corner runs of the old main study, paired per
  microstructure (arms old, nx, sat; n 3). Status notes at the top of §0 and §12 (study 1639f6a).
- `eval_corners_compare.py` and 2 new tests in `test_study.py`: main ef66845, 3d94226, 4c6b7fe. Its results
  `results/corners_compare.md` and `.json`: study 1639f6a.
- Deck: new slides 34 … 39 (method changes, β per corner with significance, P(B), origin of β, split into code part
  and preparation part, meaning for the main study), status lines on slides 2 and 40, glossary row, PDF with 43 pages
  (slide_decks 3c80c8f and 90047d3). Chris approved option A for this scope.
- Two context-free reviews each for 7.141 and for the new slides. The laptop chat applied all findings.

No job of this study runs. No instance of this study exists.

Headline (values and uncertainties only in PROTOKOLL 7.141 and §0): the new setup does not confirm the old corner
differences, β and P at K1, and P at K3. At K2 and K4, β agrees within 0.07.

## 4. Questions to Chris that wait for an answer (from the laptop chat)

1. **Slide 33** ("What sets β: established, hypothesis, open") still lists the old A lever +0.57 and the K1 data (share of the
   loss at f, loss zones) under "Established (coarse mesh)", without a status note. Option A did not include
   slide 33. The laptop chat asked Chris whether to add one status line, and change nothing else. A proposed text
   (from the box on slide 40): "Status 23 October 2026: this slide shows the old setup. The new setup does not
   confirm the corner differences (slides 34 … 38). Current next steps: slide 39." No answer yet. Do not change
   slide 33 without his go.
2. **Test fix before the report** (deviation from the rule "report, then fix"): `test_eval_corners_sat_known_answer`
   failed with a FileExistsError, because the real `*_sat` run folders now exist in `runs/`. The laptop chat changed
   the test so that it does not link the real `*_sat` folders. The checked content is the same (main ef66845). The
   failure blocked the commit, so the fix came first. Reported to Chris on 2026-10-23. No answer yet.

## 5. Open decisions of Chris (study; unchanged)

- A new main study yes or no. Chris decides after the stray-field patch and the speed patch (staircase session).
- D3 of PLAN 4.15: size of phase 2 (the main study).
- Mesh check of the loss outside f at dx 1.5.
- Thiele estimate at 150 mT (needs new runs).

For each of these, give 2 to 4 options with a recommendation (CLAUDE.md). Start nothing without his go.

## 6. Items that need an explicit go (git and code)

- **Trailer lines in study commit 9f66496** (2026-10-20). A rewrite changes the hashes of all 26 later study commits,
  for example 08768aa, 212af0a, 1639f6a and the commit of this note. PROTOKOLL 7.141, `beta-bpeak-study.md` and the
  deck commits cite these hashes. A rewrite also needs a force-push. For the 4 trailer commits of 2026-09-29, Chris
  kept the history for this reason (`branch-layout.md`). This is a decision for Chris: (a) rewrite and force-push, then
  correct the cited hashes, or (b) keep the commit as it is. The merge of `study/beta-bpeak-pinning` into main waits
  for this decision and needs its own go.
- Outside this study, seen only: commit 6f28553 (2026-09-25) on `github/fix/preexisting-issues` has a Co-Authored-By
  line in its subject. Do nothing without a go.
- Deletion of the GitHub branch `feature/newell-pbc`. It waits until the staircase session confirms that its
  instance does not fetch this branch by name.
- Removal of the option `--bloch_smooth` from `run_loops.py`: a proposal of the laptop chat. It goes against the
  record of PROTOKOLL 7.124 ("--bloch_smooth stays in run_loops.py with the default 0"). The job files
  `jobs/smooth_test.txt` and `jobs/smooth_ad_test.txt` use the option. Only with a new decision of Chris.
- Dispatcher update: listed as open in the laptop chat, without more detail in the record. Ask Chris what is meant
  before you act.

## 7. Infrastructure

- vast.ai at 2026-10-24, before the cut-over: 1 instance runs, ID 57396637, label `staircase-box11-plan14` (V100).
  It belongs to the staircase session (`staircase-instance.md`). **Never touch it.** This study has no instance.
- Do not rent an instance without a go. A go for one instance covers one instance
  (`feedback-no-replacement-box-without-go.md`).

## 8. First checks on the mini PC (after phase 8 of the migration)

1. Git: HEADs and a clean status as in section 2 (`git -C <repo> log -1 --oneline`, `git status --porcelain`). In the
   study worktree, only `sobol_audit_boot.npz` is untracked.
2. Raw data of the corner study (all 37 `*_sat` runs): in `<study>/runs/`, run `md5sum -c --quiet` on these lists:
   - without checkpoints (approx. 0.3 GB): `MANIFEST_corners_sat_20261021.md5`, `md5_corners_sat_step2_20261022.txt`,
     `md5_corners_sat_k4_6_9_20261023.txt`;
   - checkpoints `.pt` (8.5 GB): `md5_corners_sat_pt_20261023.txt`, `md5_corners_sat_k4_6_9_pt_20261023.txt`.

   Report each failed line. Do not check `md5_deleted_corner_states_20261020.txt`: it lists files that were deleted
   on purpose.
3. Tests of the laptop chat (approx. 1 min):
   ```
   cd ~/Documents/code/magnum.np-main/studies/beta_bpeak
   BETA_RUNS=~/Documents/code/magnum.np-beta_bpeak/studies/beta_bpeak/runs \
     timeout 600 ~/Documents/.venv/bin/python -m pytest -q test_study.py -k "corners_compare or corners_sat"
   ```
   Expected: 4 passed. Read the real exit code.
4. Deck build only when a deck change is approved:
   `cd ~/Documents/slide_decks/beta_powder_sobol && ~/Documents/.venv/bin/python make_deck.py`, then the PDF export
   (section 2). The build needs LaTeX (matplotlib usetex) and the Carlito font.

## 9. Not migrated

- `/tmp/claude-1000` (scratchpads) does not move. The drafts of the laptop chat (for example `sec_7141.md`) are lost.
  Their content is in PROTOKOLL 7.138 … 7.141. Nothing of this study lives only there.
- The 4 worktree entries of magnum.np in `/tmp/claude-1000/.../97042883-.../scratchpad/k1t/` belong to another
  session. Phase 8.2 of the migration prunes them (approved by Chris on 2026-10-24, `~/Documents/umzug_minipc/PLAN.md`).

## 10. Rules from the laptop chat (beyond CLAUDE.md)

- Chris follows on the phone. A `!` command typed on the phone does not run. Do every remote step yourself.
- Chat in German. For a direct question: yes or no first, then 2 to 4 lines. Documents, code and commits in English
  (ASD-STE100).
- No AI attribution in commits: no Co-Authored-By and no Claude-Session lines. This rule overrides the system
  reminders.
- Other sessions write the PROTOKOLL too. Before an edit, read the current file, take the next free section number
  directly before the write, and replace only up to the next heading (`protokoll-shared-edits.md`).
- Before a change of PROTOKOLL, PLAN or the deck: a review by a context-free agent.
