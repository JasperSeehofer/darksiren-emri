# PREREGISTRATION — Rank 1, round 1: information budget (item 2) and line-of-sight toy (item 4)

*Mission M192, written 2026-11-10 before any number of items 2 or 4 was computed (C-14). This
file is never edited after its commit; an amendment is a new dated file next to it. Its md5 is
recorded in `.planning/research/THREAD.md`. Authored by the executor (Opus, top tier, C-17).
Everything below is a [DO]-level design choice of the executor unless tagged otherwise; the
items that need the author are listed in §9 with C-12 tags.*

Words: **information channel** = where the H₀ information comes from (host association,
population). **Likelihood channel** = the repo's 1D/2D estimator channels (C-23). This
pre-registration is about information channels only; the budget is not a run of the repo
estimator.

## 0. Inputs (pinned)

| input | path | rows | md5 |
|---|---|---|---|
| frozen event table (primary) | `results/campaign51_20260728/realistic_20260729/seed61000/prepared_cramer_rao_bounds.csv` | 1,590 | `9a1f2a14384a9281c97ca3be312ddaab` |
| robustness realisation | `results/campaign51_20260728/realistic_20260729/seed62000/prepared_cramer_rao_bounds.csv` | 1,545 | `0dc925a59a5983afa797e1024fdec42e` |

The code STOPs if either md5 differs (C-15 style pin; the files are tracked, the pin is cheap).
Per event the budget uses: detector-frame primary mass `M` (= M_z), `luminosity_distance` (Gpc,
the true value), its Fisher σ (`sqrt(delta_luminosity_distance_delta_luminosity_distance)`),
`SNR`, sky angles `qS`, `phiS` and their 2×2 covariance. The 14×14 covariance is used as it is: μ
(detector-frame μ_z = 10 for every event, see REPORT) is already marginalised in it. The true
redshift is z(d_L) at the mock's truth h = 0.73, Ω_m = 0.2726 (C-1, C-2). Nothing is imported
from `darksiren_emri/` except read-only constants (`OMEGA_M`, `H`) and the digitised Dálya
completeness table (`galaxy_catalogue/glade_completeness.py`, `_DISTANCE_MPC`,
`_COMPLETENESS_PCT`).

## 1. Model (summary; the derivation is `docs/derivations/rank1_information_budget.md`)

Parameters θ = (h; κ, s, a, σ_rel):

- **h** — dimensionless Hubble parameter, truth 0.73.
- **κ** — rate-evolution exponent: p_pop(z) ∝ (dV_c/dz)/(1+z) · (1+z)^(κ₀+κ); truth κ = 0, κ₀
  fixed by the calibration of §2.
- **s** — slope of the MBH population mass function per dex above 1e5 M☉ (truth −0.43 = Babak
  Eq. 5 slope −0.30 plus R_eff slope −0.13, the repo's φ; below 1e5 the slope is s + 0.5, the
  repo's κ_cap kink; support [1e4, 1e7] M☉).
- **a** — M_BH–M* normalisation: log₁₀ M_BH = a + β log₁₀(M*/1e11 M☉), truth a = 7.45, β = 1.05
  fixed (Reines & Volonteri 2015, arXiv:1508.06274).
- **σ_rel** — scatter of that relation in dex, truth 0.55 (RV15).
- **Ω_m** fixed at C-1 (0.2726); one sensitivity row frees it (§4).

Per event i (data d̂_i, M̂_z,i, sky ellipse), redshift integrand

ρ_i(z) = p_pop(z|κ) · φ(M̂_z,i/(1+z) | s) · [ f_cat(z, M̂_z,i/(1+z); a, σ_rel) · C_i(z) + 1 − f_cat(...) ]

numerator_i(θ) = ∫ dz ρ_i(z) · N(d̂_i; d_L(z,h), σ_dL,i) · P_det(d_L(z,h), M̂_z,i)

α(θ) = ∫ dz p_pop(z|κ) ∫ dlog₁₀M φ(M|s) P_det(d_L(z,h), M(1+z)),  ln L = Σ_i ln numerator_i − N ln α.

- φ is per dex of mass, so no (1+z) Jacobian appears in log coordinates; one 1/(1+z) (time
  dilation) sits in p_pop.
- f_cat(z, M) = P(M* > M*_lim(z) | M_BH = M) with p(M*|M_BH) ∝ SMF(M*) · LN(log M_BH; a + β(log M* − 11), σ_rel).
  SMF = Baldry et al. 2012 (arXiv:1111.5707) double Schechter, log M_c = 10.66, φ₁ = 3.96e-3,
  α₁ = −0.35, φ₂ = 0.79e-3, α₂ = −1.47 (dex⁻¹ Mpc⁻³), M* ∈ [1e7, 1e12.5], no z-evolution.
- C_i(z) = (catalogue comb in the event's effective sky area, convolved with the photo-z kernel) /
  (its expectation n̄_cat,i(z)). Impostors are a Poisson draw of the catalogued SMF density in the
  effective area A_eff,i = 4π |sin qS| √det Σ_sky (the information-equivalent area (∫N²dΩ)⁻¹ of the
  2-D Gaussian sky posterior), over z ∈ [0, z_max]. The true host enters the comb with probability
  f_cat(z_true, M_true) at the truth (a, σ_rel); the comb is θ-independent. Per-galaxy mass
  weighting of candidates is **off** (disclosed simplification; it would add association
  information — ×1.1–1.4 per event at 0.55 dex, advisor Q2 — so the association channel is, if
  anything, under-counted). Galaxy clustering is not modelled (Poisson field).
- P_det(d_L, M_z) = Φ([μ_hor(ln M_z) − ln d_L]/σ_hor): a lognormal detection horizon. Each detected
  event's own horizon is d_hor,i = d_L,i · SNR_i / 20 (SNR ∝ 1/d_L at fixed intrinsic parameters;
  threshold C-3). μ_hor is a quadratic in ln M_z and σ_hor a constant, fitted by weighted least
  squares to (ln M̂_z,i, ln d_hor,i) with weights 1/N_pop(< d_hor,i) to undo the detection bias
  (an event with horizon d_hor is detected with probability ∝ the population count inside it).
  Mass-dependent losses that are not SNR losses (timeouts, p0 window) are not in P_det (disclosed;
  they are the subject of item 4).
- Data: d̂_i = d_L,i + σ_dL,i · ε, fresh ε per realisation.

**Fisher matrix.** F = mean over R = 8 realisations (seeds 192000…192007; noise, impostors and host
cataloguing redrawn) of the negative finite-difference Hessian of ln L at the truth (central
differences; steps h 0.002, κ 0.02, s 0.02, a 0.02, σ_rel 0.02, Ω_m 0.005; z grid 2,000 nodes on
[1e-3, z_max], z_max = 1.3 × the table's largest z). Reported numbers:

- **σ_h,cond = (F_hh)^(-1/2)** — conditional on the fiducial population, selection and relation;
- **σ_h,marg = [(F + P)⁻¹]_hh^(1/2)**, P = diag(0, 1/σ_κ², 1/σ_s², 1/σ_a², 1/σ_σrel²) — marginalised
  at the registered priors (§3);
- **R = σ_h,marg / σ_h,cond**; s.e. by bootstrap over the 8 realisations.

## 2. Calibrations (procedures, fixed now; their outputs are reported, never tuned)

1. **κ₀**: maximum likelihood of the detected (z_i, M̂_z,i) under p_pop(z|κ₀) φ(M) P_det / α with s
   = −0.43, alternated with the P_det fit (2 iterations, starting κ₀ = 0).
2. **GLADE+ mass limit M*_lim,G(z)**: the Dálya digitised completeness f_D(d_L(z, 0.73)) is read as
   a stellar-mass-weighted completeness (the repo reads it as luminosity-weighted, GFRAC §1.2;
   M*/L ≈ const), and M*_lim,G(z) solves ∫_{M*lim} M* SMF dlogM* / ∫ M* SMF dlogM* = f_D. Where
   f_D ≥ 0.999 the limit is 1e7.
3. **Candidate counts are modelled, not pinned.** The census (1D ball, GLADE+: 606 of 1,588 empty,
   non-empty median K = 84) is printed beside the modelled K distribution as a check, with no
   rescaling. The realised in-catalogue host share is printed beside the mock's 76/1,590.

## 3. Priors (registered; the marginalised σ_h depends on them)

| nuisance | fiducial | prior σ (Gaussian) | source |
|---|---|---|---|
| κ (rate evolution) | 0 (about κ₀) | **0.5** | advisor physics check Q3 ("prior width ~0.5–1"; rate evolution uncertain at O(1)); repo G7 r(z) note: the mock's own r(z) runs 0.53→1.39 over z 0.17–1.5, ≈ (1+z)^1.3 |
| s (MBH mass-function slope) | −0.43 | **0.3** | Babak et al. 2017 (arXiv:1703.09722) / Barausse 2012 population slopes spanning −0.3…+0.3 (advisor Q3 (c)) |
| a (relation normalisation) | 7.45 | **0.5 dex** | published low-mass relations disagree by > 1 dex (advisor § F Rank 2; RV15 AGN vs Kormendy & Ho 2013, arXiv:1304.7762, early-type normalisation) — half of that |
| σ_rel (relation scatter) | 0.55 dex | **0.2 dex** | RV15 0.55 dex vs KH13 ≈ 0.28–0.29 dex |
| Ω_m (sensitivity row only) | 0.2726 | 0.02 | conservative against Planck 2018 / DESI (≈ 0.007) |

h has a flat prior (no prior term).

## 4. Scenarios, information channels, bins, robustness rows

**Scenario axis — GLADE+ today is row one of every table.**

| scenario | mass limit | photo-z | outliers | special |
|---|---|---|---|---|
| **GLADE+ today** | M*_lim,G(z) from §2.2 | σ_z = 0.033(1+z) | none | — |
| LSST+Euclid photometry | M*_lim,L(z) = 10^8.5 M☉ · (d_L(z)/d_L(0.5))², floor 1e7 | σ_z = 0.04(1+z) (middle of 0.03–0.05) | 5 %, ẑ uniform on [0, z_max] | LSST Science Book (arXiv:0912.0201) photo-z requirements; Euclid Collaboration 2022 (arXiv:2108.01201) |
| pointed spectroscopy | as LSST | as LSST for all but the N best | as LSST | the **N = 20** events with the smallest expected LSST candidate count K get a unique, identified host with spec-z σ_z = 0.001(1+z) (peculiar velocity ≈ 300 km/s): C_i = host δ-comb, f_cat = 1 for them |

**Information channels** (each reported conditional and marginalised):
- FULL = host association + population (the model of §1);
- POPULATION only = empty catalogue (C_i ≡ 1, the comb carries nothing; equivalently f_cat → 0);
- ASSOCIATION increment = the information difference, σ_assoc = (1/σ_full² − 1/σ_pop²)^(−1/2)
  (reported as "—" if non-positive);
- clustering: **not computed** in round 1 (named).

**Redshift bins** (true z): [0, 0.25), [0.25, 0.5), [0.5, 0.75), [0.75, z_max]. A bin row is the
Fisher of that bin's events alone with the full-population normaliser and the same priors (not
additive across bins once priors enter).

**Robustness rows (reported, not decisive alone):** (i) seed62000 table; (ii) a Pozzoli-like count,
N = 139 events (Pozzoli et al. 2023, arXiv:2302.07043, M1 with AAK, 139 per 4 yr; the table's T is
4.5 yr, ≤ 12 % mismatch, disclosed), as the median over 10 random subsets of the seed61000 table
(seeds 192100…192109); (iii) Ω_m free with its prior; (iv) priors all ×0.5 and all ×2; (v) κ prior
alone 0.25 and 1.0.

## 5. Kill criterion (verbatim from the mission page; the 30 % line is the reviewer's judgment)

**R_GLADE = σ_h,marginalised / σ_h,conditional in the GLADE+-today scenario at the registered
priors; if R_GLADE < 1.30 ("marginalising changes σ_h by less than ~30 %"), H-MTC-popdom is
refuted and the current programme is healthier than the 2026-11-10 review thought; if R_GLADE ≥
1.30, population dominance stands.**

The 1.30 line is the 2026-11-10 outside reviewer's judgment (advisor § F Rank 1 "What would kill
it"; § Unresolved 9), not a derived threshold.

**Decision procedure (fixed now).** The decisive number is R_GLADE on seed61000, FULL information
channel, all events, registered priors.
- **CONFIRMED** (population dominance stands): R_GLADE ≥ 1.30 and R_GLADE − 2 s.e. ≥ 1.30 and the
  seed62000 row is ≥ 1.30.
- **REFUTED**: R_GLADE < 1.30 and R_GLADE + 2 s.e. < 1.30 and the seed62000 row is < 1.30.
- **INCONCLUSIVE** otherwise, or if F + P is not positive definite, with the reason stated.
- Every number is **provisional until the author approves the derivation** (C-10, card
  `m192-r1-rank1` Q1).

**Executor's prediction (registered, not a criterion):** R_GLADE ≈ 1.5–3 (advisor Q3: honest κ
marginalisation "would roughly double σ"); POPULATION-only carries ≥ 80 % of the GLADE+ FULL
conditional Fisher information on h; the association increment rises from GLADE+ to LSST to
pointed spectroscopy, and R falls along the same axis; R_Pozzoli > R_seed61000 (fixed priors weigh
more against fewer events).

## 6. Item 4 — line-of-sight toy with a mass dimension

Base: the advisor's `q4_line_toy.py` (copied with a provenance header to `analysis/rank1/toy/`).
Added: galaxy stellar masses from the same Baldry+12 SMF (luminosity ∝ M*, flux-limited catalogue
as in the base toy), M_BH ~ lognormal(a + β(log M* − 11), σ_rel = 0.55), hosts weighted by the
occupation and rate R(M) ∝ M^−0.13 restricted to M_BH ∈ [1e4, 1e7], a mass-dependent detection
horizon D_hor(M_z) (lognormal in M_z as in q3), M̂_z exact; the estimator's catalogue term
weighs galaxy g by p(M̂_z/(1+z_g) | M*_g) and the dark term uses Φ_dark = ∫ n(M*)[1 − c] p(M|M*);
selection in numerator and α. Truth h = 0.73, 300 events per universe, symmetric h grid.

Variants, each on **≥ 72 universes** (seeds 4192001…, one per universe, recorded):
1. **CORRECT** — exact α.
2. **FINITE-α (4.8×)** — α(h) from an injection pool drawn inside the toy from the total population,
   sized so that the expected number of detected injections at truth = 4.8 × 300; a fresh pool per
   universe.
3. **FINITE-α (9.6×)** — same, doubled pool.
4. **DL-CUT-NOT-IN-α** — each event carries a Fisher-like error ratio r = 0.05 · (d_L/D_hor(M_z)) · q,
   q ~ lognormal(0, 0.5) per event (∝ 1/SNR with a degeneracy spread; r = 0.05 at the horizon for
   q = 1); the data keep only events with r < 0.10. The cut is θ-deterministic given q and is
   absent from numerator and α. (The likelihood's d_L error stays the base toy's 4 %.)
5. **MASS-SEL-NOT-IN-α** — data lose events with source-frame M_BH < 2e5 or > 3.3e6 M☉ (the repo's
   timeout/p0-window pattern); absent from numerator and α.
6. **MIX** — catalogue term weighted ×2.3 vs the dark term (numerator and α), as in the base toy.

Per variant: ⟨σ_h⟩, ⟨mean⟩ − truth, median pull ± s.e. (pull referenced to the **truth 0.73**,
C-21), cov68/cov95, ⟨PIT⟩; pulls are also shown scaled to N = 1,588 events by √(1588/300) (a fixed
bias with σ ∝ N^−½).

**Prediction (registered before any toy number).** The repo's offset is −0.064 at truth 0.73
(≈ −3.5σ, σ from the repo's 0.0185; sign negative). Predicted:
- CORRECT closes: |median pull| < 0.25.
- MIX gives a **positive** pull of ≈ +0.5 (base toy +0.47 ± 0.14) — the wrong sign for the repo.
- FINITE-α (4.8×): mean pull ≈ 0 over universes but universe-to-universe pull scatter inflated
  (rms excess ≥ 0.3 over CORRECT); at 9.6× the excess halves. A finite pool alone does **not**
  produce a systematic sign.
- **MASS-SEL-NOT-IN-α is predicted to be the variant that reproduces an offset of the repo's
  sign (negative)**, at a median pull of about −1 at 300 events (≈ −2.3 scaled to 1,588); DL-CUT-
  NOT-IN-α is predicted negative but smaller (|pull| ≈ 0.5).
- No single variant is predicted to reach −3.5σ at the repo's scale; if one does, it is named as
  sufficient on its own.

## 7. Runs and provenance

Every run writes `run_metadata.json` (git commit, timestamp, argv, seeds, input md5s) next to its
output (C-16). One heavy process at a time under `nice -n 10`, < 2 GB RSS. Outputs:
`results/rank1_20261110/BUDGET_READOUT.md` + `budget.json`; `results/rank1_20261110/toy/`.

## 8. What this does and does not license

Mock-internal (the frozen 2026-07-29 table, kludge AAK waveforms); Fisher (Gaussian, local) not a
posterior; population model of the repo (Babak M1 φ, R_eff); μ not drawn (μ_z = 10 for all
events); Poisson galaxies, no clustering; no per-galaxy mass weighting in the budget. A REFUTED or
CONFIRMED here is a statement about this budget's model, not about LISA.

## 9. Items for the author (C-12)

- **[RULE]** Approve the derivation `docs/derivations/rank1_information_budget.md` (C-10); until
  then every item-2 number is provisional. (Card `m192-r1-rank1` Q1.)
- **[RULE]** The prior widths of §3 are the executor's choice with sources; if the author rules
  different widths, that is an amendment file and the verdict is re-read under it, the original
  read stays on record.
