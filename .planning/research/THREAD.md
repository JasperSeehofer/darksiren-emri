# Thread state — EMRI → H₀ dark-siren inference

*Layer-2 thread-state contract (vault: `wiki/analyses/research-routine-packs-spec.md`). Read by
the vault's `/cartographer` bridge. Keep current: same-day update after every experiment run
(vault: `wiki/meta/research-operating-manual.md`, "per numerical experiment / run"). Vault mirror
of the hypothesis rows: `wiki/meta/hypothesis-ledger.md`.*

> **Created 2026-09-21 by mission M65** (first research project wired into the manyfold vault).
> Vault entity page: `wiki/entities/master-thesis-code.md` — note the vault slug is
> `master-thesis-code`, deliberately decoupled from this repo's directory name
> (see `wiki/meta/darksiren-emri-rename-plan.md`). Experiment/hypothesis prefix is **MTC**, the
> prefix already in use in the vault ledger — not a new one coined here.
>
> **Continuity warning.** `STATE.md` (repo root, "Last updated: 2026-08-12") is **stale** and must
> not be read as the thread's frontier: `main` HEAD is 2026-08-26, while the live m1 work sits on
> the unmerged branch `fix/p32d-classg-venue-repair` (HEAD `d201d9e4`, 2026-09-04). The whole
> `graph1_20260901/` tree — 77 commits — **does not exist on `main` at all**. STATE.md still says
> T-1 is "still unrun", which remains true, but it predates the entire m1 exercise.

## Hard core

**[PROPOSED — awaiting author ratification as a [RULE]; not ratified by M65.]**

The closure of the dark-siren H₀ estimator at the injected truth is a consequence of
**generator-consistent normalisation** (n̂_w = W_cat/V_f + D_gen, plus z-resolved survival
S(d_L|z)) — that is, of the physics of the selection function — and **not** of any
truth-referencing constant, anchor, or tuning introduced by the estimator.

*Source of the formulation:* `H-MTC-genmarg` in `wiki/meta/hypothesis-ledger.md`, promoted to
thread level as kill criterion (1). The `project-init` skill deliberately leaves the hard core as
a fill-in stub because it is the author's to set; the statement above is a **mirror of the
author's own already-written ledger row**, offered for ratification, not a new claim invented by
the mission.

## Live hypotheses (mirror of vault hypothesis-ledger)

| id | statement | discriminating experiment | status |
|---|---|---|---|
| H-MTC-genmarg | The deep-venue closure at truth comes from generator-consistent normalisation plus z-resolved survival S(d_L\|z), not from any truth-referencing constant | **T-1** blind alternative-truth mock: seal h_inj, re-run the production stack, check the MAP tracks the sealed truth | **live** (opened 2026-07-26) — T-1 **unrun**, see Active experiment |
| H-MTC-realdata | The precision the δ-kernel host-z numerator buys survives the move from mock to real data once the photo-z kernel replaces point/point pairing | re-derive the photo-z kernel, then re-measure σ_h against the mock-internal baseline on the same venues | **live** (opened 2026-07-28) — unstarted |
| H-MTC-norm | The H₀ rail is a curable normalisation artifact, not information starvation | run the production estimator under `{global, local_ratio, volume_deconv}` | merged 2026-07-02 (named winner superseded by H-MTC-genmarg) |
| H-MTC-massch | The 2-D mass channel tightens H₀ by ~2× (4%→1.8%) | forecast with realistic host-mass scatter | pruned 2026-07-01 — information-starved |
| H-MTC-vtransfer | The ball-venue calibration DEFECT is a venue caricature, not a property of the production estimator | production-matched venue-transfer campaign | **REFUTED** 2026-08-13 by its own decider (TRANSFER-CONFIRMED, author-ratified); vault status change to `pruned` is **proposed, not applied** — a prune is a kill and the human confirms kills |

## Kill criteria

Thread-level, **written 2026-08-11 by the author** (vault: `wiki/meta/hypothesis-ledger.md`,
[[proposal-queue|PROP-20260811-01]]). Quoted, not paraphrased:

1. **Anti-tuning** — if T-1 (blind sealed-h mock) closes anywhere other than the sealed `h_inj`,
   the `generator_marginal` + `--pdet_z_resolved` adoption **reverts** and the closure claim is
   **withdrawn, not patched**.
2. **Precision** — if the real-data photo-z kernel degrades σ_h to the host-known limit (×3.3–×6.8)
   or worse, every precision claim is **restated as mock-internal** and the manuscript's headline σ
   is withdrawn.
3. **Time-bound** — if no blind-mock verdict exists by **2026-09-23**, the thread ships as a
   **methods** paper with the closure reported as *unblinded*, rather than as an H₀ measurement.
   *(Date borrowed from the workspace copy-off deadline already on the books, not derived from the
   physics — amendable.)*

**Status of (3) as of 2026-09-21:** extended by author ruling on card `overnight-report-2026-09-21`
(*"Explicit extension: run T-1 + adjudicate m1 firs[t]"*), conditioned on running T-1 **and**
adjudicating the m1 docket first. M65 has adjudicated m1 (see Current claims); **T-1 could not be
run** (see Active experiment), so the extension's stated precondition cannot be met as written.
Returned to the author as a `go-kill` review card by M65.

## Current claims

| claim | status | evidence |
|---|---|---|
| The m1 sealed-mock's registered cell disposition | **UNRESOLVED-AS-REGISTERED** — the pre-registration does not cover a posterior on which TUNED and NOT-TUNED-AT-RAIL both hold; no tie-break invented | M65 adjudication, `.planning/missions/M65-darksiren-emri-into-manyfold/REPORT-R1.md`; raw record `results/campaign51_20260728/realistic_20260729/graph1_20260901/exec/r-sealed-mock/READ_RECORD_m1.md` (branch `fix/p32d-classg-venue-repair`) |
| m1 posterior is centred at 0.628–0.631, i.e. 0.0423 (2D) / 0.0389 (1D) **below** the 0.67 truth, with 71 %/62 % of mass at or below 0.63 | verified twice (chair + M65 independent verifier); cell-free, stands regardless of the R25 ruling | `READ_RECORD_m1.md` §4; M65 `VERIFY-R1.md` §1 |
| m1 makes **no anti-tuning claim** and does not discharge T-1 | registered, explicit | `D_SEALED_REGISTER_DOSSIER.md` §2 (*"the truth is public (no anti-tuning claim can bank from (m1))"*), §1.1 (*"explicitly NOT the T-1 verdict"*) |
| Anti-tuning (kill criterion 1) | **untested** — T-1 unrun | — |
| Precision on real data (kill criterion 2) | **untested** — photo-z re-derivation unstarted | — |
| Every decisive m1 number is prose-only | **defect** — the underlying posteriors and the scorer are not in version control | M65 `VERIFY-R1.md` §3 T1 |

## Active experiment

**T-1 (blind sealed-h mock) — BLOCKED BY TWO INDEPENDENT GATES.** Not running, not scheduled, and
*not* substituted for.

In this campaign's framing T-1 is the **sealed (m2) run**. It has never been launched: the only
sealed-mock exec directory on the branch is `r-sealed-mock` (m1); no `(m2)` directory exists.

**Gate 1 — no cluster credential (infrastructure).** T-1 requires a fresh GPU EMRI simulation
campaign at a newly sealed `h_inj` on bwUniCluster 3.0 (SLURM, `gres=gpu:1`, CUDA 12 —
`cluster/simulate.sbatch`), then a CPU eval array
(`results/redteam_20260726/PHYSICS_METHODOLOGY_REVIEW.md` §4: *"Cost: one simulation campaign +
one eval array"*). The brain host holds no SSH key and no university VPN; the mission built to
solve that, `M38-cluster-session-unlock`, is itself `stage: idea` and unbuilt. The redteam review
rules the CPU-only "cheap partial" (rescaling an existing CSV's d_L) **"NOT a substitute"** — it
*"confirms the §1 mechanism; it does not test the selection model."*

**Gate 2 — the author's own launch bar (governance).** Docket **R27** states: *"no (m2) launch
until ratified"*, pending an erratum that adds a both-fire precedence rule and a stated treatment
for σ_measured < σ_anchor. **This gate is independent of the credential**: even with cluster access
tonight, T-1 could not legitimately launch. Ratifying R27 is an author action that unblocks T-1's
governance side regardless of when the credential lands.

**Do not book m1 as T-1.** m1 is a sealed *mock* at an alternative truth, sharing T-1's design but
not its blinding: its truth (h = 0.67) is public and legible in the run directory name
(`run_20260729_seed64000_h0p67`). The register says so in its own words, repeatedly.

*No `EXP-MTC-*` row is pre-registered in `wiki/meta/agentic-experiments-research.md` yet; the next
experiment to pre-register is T-1 itself, once the credential lands.*

## Claim history

`CLAIMS.jsonl` (commission-compatible) is created by the first `/commission --research` run on
this thread; until then the table above is the claim record. Note the repo also maintains a
long-running append-only record at
`results/campaign51_20260728/realistic_20260729/gate_b_20260730/BIAS_HISTORY_LEDGER.md`
(row #387 is the m1 read), which predates this file and is not superseded by it.

## Last verification

**2026-09-21 — M65 independent verification pass (`VERIFY-R1.md`).** Fresh-context adversarial
verifier, dispatched before any adjudication draft existed. Verdict: the four m1 numbers
(2.73σ / 2.58σ / 4.39 % / 2.59 %) **all reproduce exactly** against the record; the
seal-before-result ordering is **evidenced by git** (threshold block byte-identical across all
four dossier commits, hash `c908790d…`, frozen ~11 h before the first array task); but the
"both cells fire" reading is **CONTESTED, not a fact**, and the underlying posterior data is
**absent from version control**.
