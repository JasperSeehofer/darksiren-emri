# M192 — REPORT, round 1 (executor)

**R_GLADE verdict (one line):** **INCONCLUSIVE** under the registered rule, because F + P is not
positive definite along (h, κ) on seed61000 and seed62000. **Provisional** until the C-10
approval. The post-hoc OPG diagnostic, not registered, gives R ≈ 2.7 against the 1.30 line.

## Least-confident decisions (read these first)

1. **The registered Fisher estimator (observed information at the truth, on the real table) is
   not positive definite in the GLADE+/LSST scenarios.** The model's p(z | κ) and fitted P_det
   do not reproduce the table's z > 0.5 tail (κ₀ = −1.35, quadratic horizon). The observed
   curvature in (h, κ) is a saddle there. I applied the registered rule (INCONCLUSIVE) and
   added two labelled post-hoc diagnostics: the outer product of scores, and the eigen-structure.
   An expected Fisher matrix on model-drawn data would have been the better registration. That is
   my design error, now round-2 item 2.
2. **The GLADE+ mass limit makes catalogued hosts vanish (0.01 %).** I read the Dálya completeness
   as mass-weighted and held it flat beyond 800 Mpc, as the repo's `glade_completeness.py` does.
   Beyond z ≈ 0.17 the limit is then M* ≈ 1e11. That follows the registered procedure and is
   physically defensible (advisor Q2: ~2 % at a 1e10 limit). It is far below the mock's 4.8 %,
   which is a mass-independent Bernoulli draw. The GLADE+ association channel is therefore ~0 by
   construction of this reading.
3. **Toy FINITE-α.** The registered variant (whole α from a population-drawn pool) is biased by
   −0.129 at both pool sizes. The INF-POOL / FINITE-DARK split that attributes this to "population
   α vs realised-catalogue numerator" was added after a development run on non-registered seeds
   (900–907). The diagnostics are labelled post hoc. The registered rows are reported unchanged.
4. **Numerical fix after a first look.** The first budget run used a per-event z window that
   under-resolved events with σ_dL/d_L ≈ 1e-3 (one event at −2.3e4 curvature). I stopped it,
   widened the window and re-ran. The registered "2,000-node global grid" was never implemented
   as such (per-event 161-node grids instead). Both deviations are in the readout's deviations
   section. No number had been read against the criterion: the first run's R_GLADE was undefined too.
5. **Code at run time vs committed code.** The budget and toy runs used the working tree.
   Both `run_metadata.json` files record HEAD `14f7fd23` (item 1) with the code uncommitted. The committed files
   differ from the run-time code only by `ruff format` and three `np.asarray` wrappers for mypy.
6. **Base test check (d).** The first full run hung for > 30 min on a futex-deadlocked forked
   worker. `pkill` was refused by the auto-mode classifier, so I stopped my own task with
   TaskStop. A re-run with a log passed: **2,033 passed, 31 skipped, 0 failed** (6 min 9 s).
7. **Card Q2 (round 2)** offers three options. The suggested answer is "expected Fisher + joint
   likelihood", because it turns the INCONCLUSIVE into a number. The note says the numbers alone
   favour Rank 3 (0.5 % from 20 spectroscopic hosts, R ≈ 1).

## Budget

The executor cannot see the meter. Work done: one Opus orchestrator session that wrote the
derivation, pre-registration, budget code, toy extension, note and cards itself, plus 3 Sonnet
subagents (FEW frame read, paper edits, card-rule digest). No top-tier subagents were used. Runs:
toy 72 universes (≈ 13 min), budget twice (≈ 12 min each), one discarded and two development runs,
and the full test suite twice. Wall clock: started 22:25:56 CET (CLOCK.txt); REPORT at ≈ 23:50,
≈ 1 h 25 min, no MEMORY-RULE pauses (load stayed < 10). My guess is well inside the 8 M declaration,
but the overseer's meter is authoritative.

## Items

| item | status | commit |
|---|---|---|
| 0 start checks | (a) M162 ancestor ✓ (b) tree `71acb98b…` ✓ (c) seed61000 1,590 rows, md5 `9a1f2a14384a9281c97ca3be312ddaab` ✓ (seed62000 `0dc925a5…`, 1,545 rows) (d) pytest green: 2,033 passed (see decision 6) | — |
| H hard core | ratified wording verbatim in THREAD.md § Hard core, old text kept as dated superseded note | `88a7a46b` |
| P pre-registration | `PREREGISTRATION_RANK1.md`, md5 `90dc4200f95f82c8a64b26f302f739c6` in THREAD.md, committed before any item-2/4 number | `f791c86c` |
| 1 paper | M193 grep (`/usr/bin/grep`, page Done criterion 2) → **0 lines**. Four remnants cleared (results.tex N_det=223, p_det=1 statements in results.tex and discussion.tex, Δh~0.02) with `\pending{Rank-1 budget}`. appendix_parameters.tex aim reworded to ruling 3. Title and P12 → card. No TeX on host: not compiled | `14f7fd23` |
| 2 budget | derivation, ledger row (presented / PENDING-AUTHOR), `analysis/rank1/budget.py`, 8 limiting-case tests, readout + JSON. GLADE+ first in every table | `e55e5705` |
| 3 empty-catalogue control | **round 2** (named; nothing fetched) | — |
| 4 toy | verbatim copy with provenance header (byte-identical body, md5 `d3e8c73c…`), `q4_mass_toy.py`, six variants × 72 universes, readout | `ee90919f` |
| 5 note | `docs/RANK1_NOTE.md`, THREAD.md Rank-1 block, CONVENTIONS C-28…C-31 | `201c234b` |
| cards, REPORT, CLOCK | three cards staged (size check and lint pass), this report | this commit |

Freeze: `git rev-parse HEAD:darksiren_emri` = `71acb98b0fd8a0f9d3aad9595eb9fc084cae0b6d` (asserted
before this commit). `main` (`9b028075`) was merged in: already up to date. No push, no merge to main.

## Key numbers (provisional; seed61000)

- GLADE+ today, full information channel: σ_h,cond = 0.050 (6.9 %). Population-only is identical;
  the association increment is negligible. σ_h,marg is undefined (registered). OPG R = 2.73
  (seed62000: 2.77). At N = 139, R = 6.4.
- LSST+Euclid photometry: σ_h,cond 0.049. Association increment σ = 0.29. OPG R = 2.62.
- Pointed spectroscopy (20 smallest-volume events): σ_h = 0.0037 (0.51 %), R = 1.006.
- Toy: CORRECT +0.12 ± 0.15 (median pull, ref. truth 0.73). FINITE-α 4.8× / 9.6×: −0.129 in h
  (median −2.8 / −2.7). MASS-SEL −0.39 ± 0.27. DL-CUT +0.23 ± 0.16. MIX +0.39 ± 0.20.

## μ frame convention (reading only)

FEW's `Pn5AAKWaveform` is generated with `frame="detector"`
(`darksiren_emri/waveform_generator.py:98-104`). The FEW 2.0.0rc1 docstring
(`few/waveform/waveform.py:52-55`) says detector-frame waveforms "are left alone". Masses enter
only via GM/c³ in observer time. The repo lifts M (`parameter_space.py:268`, M_z = M(1+z)) and
not μ (`:278` passes `mu.value` raw; `cosmological_model.py:180-195` sets μ = 10, `is_fixed`).
μ = 10 in all 1,590 rows. **So every event has detector-frame μ_z = 10 M☉ (source-frame
10/(1+z)), and μ_z carries no redshift information in the mock.** Recorded as CONVENTIONS C-31.

## Round-2 list

1. **Empty-catalogue control (item 3)**, the Gray et al. 2023 Fig. 6 analogue: the 0.73-anchor
   configuration, 2D primary and 1D replicate likelihood channels, an input-level empty catalogue
   (zero-row `--observed_catalogue` with a valid sidecar; a code change is R24 territory), the
   pool `mix200k_20260728` pinned by list-md5 `75f4030d…` (C-15), and a timed subsample probe
   first. Input source — emri2 before **2026-11-22**, the cluster window, or the Storage Box — is
   a [RULE] for the author.
2. Expected Fisher on model-drawn data (registered amendment file) to give R_GLADE a defined value.
   Then the joint hierarchical likelihood (card Q2).
3. A reading check of the repo's α construction against the toy mechanism: population-level α
   from injections vs a realised-catalogue numerator. Reading only, R24.
4. Paper: apply the title and P12 rulings (card `m192-r1-paper`). Fill the `\pending{Rank-1 budget}`
   markers once the derivation is approved.
5. Optional: a smooth-P_det toy variant (removes the step-function α artefact at small pools).

## Vault rows for the overseer to apply after landing (the executor never edits the garden)

- `wiki/meta/agentic-experiments-research-physics.md`: **EXP-MTC-RANK1** — Fisher information
  budget + mass-dimension line-of-sight toy; prereg md5 `90dc4200…`; result INCONCLUSIVE (registered;
  F+P not PD along h–κ), OPG diagnostic R ≈ 2.7, provisional pending C-10; toy: population-level α
  vs realised-catalogue numerator reproduces the offset's sign and size.
- [[hypothesis-ledger]] **H-MTC-popdom**: status stays **live**; evidence row "R_GLADE INCONCLUSIVE
  (registered), diagnostics indicate R ≈ 2.7 > 1.30; decisive read in round 2 (expected Fisher)".
- [[research-health]]: row 2026-11-10 "darksiren Rank 1 round 1: pre-registered, budget INCONCLUSIVE
  by its own rule (estimator not PD), toy localises a normaliser mechanism; C-10 approval pending".

## Items for the author (C-12)

- [RULE] Approve the derivation (C-10) — card `m192-r1-rank1` Q1.
- [RULE] Round-2 content — card `m192-r1-rank1` Q2.
- [RULE] Title; P12 scatter wording — card `m192-r1-paper` Q1, Q2.
- [RULE] Source of the empty-catalogue control inputs (round 2) — next round's card.
- [DO] Land — card `m192-land`.

Rubric (I1–I8, A1–A4) not run by the executor: it needs the parallax app, and this mission writes
nothing there. Not claimed.
