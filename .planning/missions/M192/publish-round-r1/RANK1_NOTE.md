# What LISA EMRI dark sirens measure — Rank-1 note, round 1 (seed of the methods note)

*Mission M192, 2026-11-10. **Provisional**: the derivation behind every budget number awaits the
author's approval (C-10, card `m192-r1-rank1` Q1). Mock-internal, Fisher-level. Pre-registration
`.planning/missions/M192/PREREGISTRATION_RANK1.md` (md5 `90dc4200…`), committed before any number.*

## The question (ruling 3)

How much H₀ information do LISA EMRI dark sirens carry once the population, the selection and the
black-hole–galaxy relation are marginalised? Where does it come from — host association or the
population (information channels)? How does it depend on the catalogue — GLADE+ today, LSST+Euclid
photometry, pointed spectroscopy (ruling 4)? H₀ is one output; the reference forecasts assume
catalogues that may not hold the typical host, and this note quantifies what remains.

## Budget (frozen seed61000 table, 1,590 events; σ_h in units of h, truth 0.73)

Full numbers, bins and robustness rows: `results/rank1_20261110/BUDGET_READOUT.md`.

| scenario | information channel | σ_h conditional | σ_h marginalised (registered) | R (registered) | R, OPG diagnostic (post hoc) |
|---|---|---|---|---|---|
| **GLADE+ today** | full | 0.050 (6.9 %) | undefined (F+P not pos. def.) | — | **2.73** |
| GLADE+ today | population only | 0.050 (6.9 %) | undefined | — | 2.73 |
| GLADE+ today | association increment | 1.16 (negligible) | — | — | — |
| LSST+Euclid photometry | full | 0.049 (6.8 %) | undefined | — | 2.62 |
| LSST+Euclid photometry | association increment | 0.29 | — | — | — |
| pointed spectroscopy (20 events) | full | **0.0037 (0.51 %)** | 0.0037 | 1.006 | 1.004 |

seed62000 (1,545 events) agrees: GLADE+ σ_h,cond 0.055, OPG R 2.77, spectroscopy 0.0031 (R 1.005).
At a Pozzoli-like count (N = 139), the registered matrix is positive definite. There R_GLADE = 6.4,
R_LSST = 3.9 and R_spec = 1.004 (σ_h,cond 0.165 for GLADE+, 0.014 for spectroscopy). Freeing Ω_m
(prior 0.02) changes σ_h,cond by < 3 %.

**Where the information is.** With GLADE+ the host is catalogued for ≈ 0.01 % of events once
cataloguing depends on the host's stellar mass. The mock's 4.8 % comes from a mass-independent
Bernoulli draw. Essentially all the conditional information is population information: the d_L
and M_z distribution read against an assumed rate evolution and mass function. LSST depth catalogues
≈ 78 % of hosts, but ~3×10⁴ candidates per box at σ_z = 0.04(1+z) leave a small association
increment. Spectroscopic redshifts of the 20 smallest-volume events give 0.5 % on their own, and
they do not care about the population (R ≈ 1).

## Verdict against the kill criterion (registered before any number)

R_GLADE = σ_h,marg/σ_h,cond in the GLADE+-today scenario at the registered priors; R < 1.30 refutes
H-MTC-popdom.

**INCONCLUSIVE (registered), provisional.** In the GLADE+ and LSST scenarios the observed
information matrix plus prior precision has a negative eigenvalue along the (h, κ) direction
(rate evolution), on both tables. The model's redshift distribution and selection do not reproduce
the table's z > 0.5 tail, and the per-bin observed information is negative there. So σ_h,marg is
undefined, and the registration names this case INCONCLUSIVE. Three things point the same way, and
none of them is registered:
- the outer-product-of-scores matrix (PSD by construction) gives R ≈ 2.7 on both tables;
- the registered matrix at N = 139 gives R = 6.4;
- the 0.25 ≤ z < 0.5 bin gives 1.33 / 1.98.

All lie above 1.30. Population dominance is therefore *indicated, not established*. A registered
number needs an expected Fisher matrix computed on data drawn from the model itself (round 2).

## What the toy says about the ~3.5σ offset

`results/rank1_20261110/toy/TOY_READOUT.md` (72 universes × 300 events, six registered variants):
- The correct likelihood with a mass dimension closes at a 2 % in-catalogue share (median pull +0.12 ± 0.15).
- **A normaliser α taken from a population-drawn injection pool** gives a large negative offset at
  4.8× and 9.6× alike: ⟨mean⟩ − truth = −0.129, median pull ≈ −2.7, ≈ −6 scaled to 1,588 events.
  The numerator's catalogue term uses the realised catalogue, so this α is mismatched to it. The
  post-hoc split shows an infinite pool is equally biased (−0.136). Pool size alone matters much
  less: with only the dark part of α from the pool, the offset is −0.31σ at 4.8× and −0.01σ at
  9.6×. This is the repo's sign and size, and the advisor's ranking #1 (leg-inconsistent
  normalisation) made concrete.
- Mass-dependent losses absent from α (timeouts/p0 window) give the right sign but only −0.4σ (−0.9 scaled).
- The σ_dL/d_L cut and the ×2.3 mixture weight give *positive* pulls; neither can explain the offset.
- Against the registration: CORRECT and MIX met; MASS-SEL right sign, 40 % of the predicted size;
  DL-CUT and "a finite pool has no systematic sign" not met.
- Whether the repo's α has this mismatched structure is a round-2 check on the code (reading, not changing — R24).

## Round 2 (named)

1. **Empty-catalogue control** (Gray et al. 2023 Fig. 6 analogue): the existing estimator on the
   seed61000 events with an empty catalogue, same configuration as the 0.73 anchor, 2D primary and
   1D replicate likelihood channels. It needs the pinned pool `mix200k_20260728` (list-md5
   `75f4030d…`). Its source — emri2 before 2026-11-22, the cluster window, or the Storage Box — is a
   round-2 decision.
2. A registered **expected** Fisher matrix (data drawn from the model; κ₀ and P_det re-fitted) to
   turn the INCONCLUSIVE into a number, then the joint hierarchical likelihood on the full stack.
3. A read of the repo's α construction against the toy's mechanism (population-level α vs realised-catalogue numerator).

## What this does and does not license

It licenses: inside this mock, Fisher-level, (i) GLADE+ adds essentially no association
information beyond the population; (ii) a handful of spectroscopic hosts carry more H₀ information
than all dark events together, and that information is population-independent; (iii) a mismatched
normaliser alone can produce the repo's offset.

It does not license: any statement about LISA data, a posterior (Fisher is local), or a
marginalised σ_h (undefined at the registered setting). It also says nothing about the μ_z
information channel: in the mock μ is not drawn. FEW's Pn5AAK waveform is detector-frame, so every
event has μ_z = 10 M☉ (source-frame 10/(1+z)), and μ_z carries no redshift information here (C-31).
Waveforms are AAK kludge waveforms. Other limits: no galaxy clustering, no per-galaxy mass
weighting, and an SMF without evolution.

## Which follow-on the numbers favour

**Rank 3 (golden EMRIs and pointed spectroscopy) is favoured.** Twenty spectroscopic hosts give
0.5 % with R ≈ 1, against ~7 % (conditional) from all 1,590 dark events. **Rank 2 (the reverse
problem) is disfavoured at GLADE+ and LSST depth.** The association signal is swamped: K ≈ 3×10⁴
candidates per box, and 0.01 % catalogued hosts with GLADE+. It revives only with spectroscopy.
**Rank 4 (μ_z)** cannot be judged on this table, which has no drawn secondary-mass spectrum.
