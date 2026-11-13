# Rank-1 information budget — derivation (C-10 presentation)

*Mission M192 scope item 2, 2026-11-10. Status: **presented / PENDING-AUTHOR** (ledger row in
`docs/gates/PHYSICS-GATE-LEDGER.md`). Every number computed from this model is provisional until
the author approves it (card `m192-r1-rank1`, Q1, [RULE]). The code lives in `analysis/rank1/`,
touches no physics-trigger file and changes nothing under `darksiren_emri/` (R24, C-26); the gate
is run because the model is new physics, not because a trigger file changes.*

The five C-10 items: (1) the formula (old → new), (2) reference, (3) dimensional analysis, (4)
limiting cases, (5) regression tests. "Old" here is: no information budget existed; the closest
precursors are the advisor's q3 toy (one-parameter dark spectral siren) and the repo's
`docs/SIGMA_Z_SIGMA_M_FORECAST.md` (synthetic forecast, pre-redesign).

## 1. The model

### 1.1 Population, data, selection

EMRI rate density in redshift and source-frame log mass:

  dN/dz dlog₁₀M ∝ p_pop(z|κ) φ(M|s),  p_pop(z|κ) = (dV_c/dz) (1+z)^(κ₀+κ−1),

the (1+z)⁻¹ being time dilation of the source-frame rate. φ is the repo's EMRI-weighted MBH mass
function per dex, φ ∝ 10^{s(log M − 5)} above 1e5 M☉ and 10^{(s+0.5)(log M − 5)} below, on
[1e4, 1e7] (Babak et al. 2017 Eq. 5 slope −0.30 times R_eff ∝ M^−0.13; truth s = −0.43).

Data per event: d̂ (Gaussian, σ_dL from the event's Fisher covariance), sky position (Gaussian,
2×2 block of the covariance), M̂_z (σ_Mz/M_z ~ 1e-7: a δ-function). The detector-frame mass is
M_z = M(1+z). In log coordinates the δ-function Jacobian is a data constant, so the source-frame
mass at trial redshift z is M̂_z/(1+z) and φ enters as φ(M̂_z/(1+z)) with no extra factor.

Selection is θ-deterministic (noise-free SNR, as in the mock), so it sits in the numerator and in
α (advisor physics check Q1; Mandel, Farr & Gair 2019, arXiv:1809.02063):

  P_det(d_L, M_z) = Φ([μ_hor(ln M_z) − ln d_L]/σ_hor).

### 1.2 Host: catalogued or dark

Host M* given the BH mass, by Bayes on the stellar mass function (SMF, Baldry et al. 2012) and the
relation log M_BH = a + β log(M*/1e11) with lognormal scatter σ_rel (Reines & Volonteri 2015):

  p(log M*| M) ∝ SMF(M*) · N(log M; a + β(log M* − 11), σ_rel).

Probability the host is catalogued at redshift z, for a mass-limited catalogue M* > M*_lim(z):

  f_cat(z, M; a, σ_rel) = ∫_{log M*_lim(z)} p(log M*|M) dlog M*.

This is the one population used by both legs: catalogue share f_cat and dark share 1 − f_cat sum to
the population (advisor Q1 coherence condition, with the population mass function φ in place of
SMF⊗relation — they agree to ~20 %, advisor Q2).

### 1.3 Per-event likelihood

The Gray-type mixture (Gray et al. 2020, arXiv:1908.06050), written as a redshift integrand:

  ρ_i(z|θ) = p_pop(z|κ) φ(M̂_z,i/(1+z)|s) [ f_cat(z, M̂_z,i/(1+z)) C_i(z) + 1 − f_cat(z, M̂_z,i/(1+z)) ]

  L_i(θ) = ∫ dz ρ_i(z|θ) N(d̂_i; d_L(z,h), σ_dL,i) P_det(d_L(z,h), M̂_z,i) / α(θ)

  α(θ) = ∫ dz p_pop(z|κ) ∫ dlog M φ(M|s) P_det(d_L(z,h), M(1+z)).

C_i(z) is the catalogue contrast in the event's sky area: Σ_g K(z; ẑ_g) / n̄_cat,i(z), with K the
photo-z kernel and n̄_cat,i(z) = n_cat(z) A_eff,i (dV_c/dz dΩ) the expected catalogued density, so
⟨C_i⟩ = 1 over impostor realisations: a field of impostors carries no information, a true host
adds a bump. The sky enters only through A_eff,i = (∫ N_sky² dΩ)⁻¹ = 4π |sin θ| √det Σ_θφ, the
information-equivalent area of the 2-D Gaussian (each catalogued galaxy in the sky ellipse
weighted equally inside A_eff). α contains no catalogue: catalogued plus dark is the whole
population.

### 1.4 Information by channel

- **Population only (empty catalogue):** C_i ≡ 1 ⇒ ρ_i = p_pop φ: the d̂ and M̂_z distribution
  against the population × selection (a spectral/"distance-distribution" siren).
- **Full:** with the comb.
- **Association increment:** I_assoc = I_full − I_pop on h (conditional or marginalised).

Fisher: F_jk = −E[∂_j ∂_k ln L], estimated as the mean over R noise/catalogue realisations of the
finite-difference Hessian at the truth. σ_h,cond = F_hh^(−1/2); σ_h,marg = [(F + P)⁻¹]_hh^(1/2) with
P the diagonal Gaussian prior precision of the nuisances.

### 1.5 Nuisance derivatives (analytic forms, used by the tests)

- ∂_κ ln ρ_i = ln(1+z); ∂_κ ln α = ⟨ln(1+z)⟩_det (detected-population mean).
- ∂_s ln ρ_i = (log₁₀ M̂_z − log₁₀(1+z) − 5) ln 10 (above 1e5; same below); ∂_s ln α = the
  detected-population mean of the same.
  For a pure power law (no kink, φ normalised) and mass-independent selection,
  ∂_s ln ρ_i = ln10·(log₁₀M̂_z,i − 5) − ln(1+z) − ∂_s ln Z(s): a per-event constant minus
  ln(1+z), while ∂_s ln α = −∂_s ln Z(s) is h-independent. The z-dependent part of the s-score is
  exactly minus the κ-score's numerator part, so the h–s and h–κ Fisher elements obey the
  identity **F_hs = −F_hκ + N ∂_h∂_κ ln α**: in the h-direction the mass slope acts as a
  rate-evolution prior up to the selection term (advisor Q3 "a rate-evolution prior in
  disguise"). The kink at 1e5 and the M_z-dependent selection break the identity.
- ∂_a, ∂_σrel enter only through f_cat: ∂ ln ρ_i = ∂f_cat (C_i − 1)/[f_cat C_i + 1 − f_cat]. For
  C_i ≡ 1 they vanish identically: an empty catalogue carries no information on the relation.
- ∂_h: through d_L(z,h) = (1+z) D_c(z)/h only, ∂ ln d_L/∂ ln h = −1 at fixed z.

## 2. References

Gray et al. 2020 (arXiv:1908.06050) completion likelihood; Mandel, Farr & Gair 2019
(arXiv:1809.02063) selection; Babak et al. 2017 (arXiv:1703.09722) EMRI population; Reines &
Volonteri 2015 (arXiv:1508.06274) relation; Baldry et al. 2012 (arXiv:1111.5707) SMF; Dálya et al.
2022 (arXiv:2110.06184) GLADE+ completeness; Barack & Cutler 2004 (gr-qc/0310125) sky-area
convention; Gray et al. 2023 (arXiv:2308.02281) empty-catalogue control (round 2).

## 3. Dimensional analysis

- h dimensionless; d_L = (1+z) D_c(z)/h with D_c in units of c/(100 km/s/Mpc) = 2997.92458 Mpc
  (h-free); σ_h carries the dimension of h (none).
- N(d̂; d_L, σ) has units Gpc⁻¹ and p_pop, φ, f_cat, C_i, P_det are densities per unit z, per dex,
  or dimensionless; numerator and α differ by one factor Gpc⁻¹ per event, a data constant that
  drops from every derivative. ln L is dimensionless; F_hh has units h⁻².
- A_eff in steradians × (dV_c/dz dΩ) in Mpc³ × n_cat in Mpc⁻³ per dex integrated over dex ⇒ K
  dimensionless (h-factors: V ∝ h⁻³ and the SMF in h=0.7 units; the comb is evaluated at the
  truth, h-independent by construction).

## 4. Limiting cases (each a test in `darksiren_emri_test/rank1/`)

1. **Complete catalogue, unique hosts, no selection** (C_i = δ at z_host, σ_z → 0, P_det ≡ 1):
   I_hh = Σ_i (d_L,i/h)²/σ_dL,i² ⇒ σ_h/h = [Σ_i (d_L,i/σ_dL,i)²]^(−1/2) → σ_dL/d_L/√N for equal
   errors.
2. **Unique host with photo-z σ_z:** σ_lnh² = (σ_dL/d_L)² + (γ σ_z/(1+z))², γ = d ln d_L/d ln(1+z)
   = 1 + (1+z)/(E(z) D_c(z)) (in c/H₀ units); γ(0.47) = 3.76 at Ω_m = 0.2726 (advisor correction
   (b)).
3. **Empty catalogue → population only:** C_i ≡ 1 gives F_full = F_pop and F_a· = F_σrel· = 0.
4. **Zero prior width → conditional:** σ_marg → σ_cond as all prior σ → 0.
5. **Slope–evolution identity:** without the kink and with mass-independent selection,
   F_hs = −F_hκ + N ∂_h∂_κ ln α (§1.5).

## 5. Regression tests

`darksiren_emri_test/rank1/test_budget.py`: one test per limiting case above, plus input pins
(md5 of the event table) and a determinism test (same seed ⇒ same F). No pre-existing numerical
result exists to regress against (new module); the limiting cases are the regression anchors.

## 6. Disclosed simplifications

Fisher (local Gaussian) not posterior; per-galaxy mass weighting off; Poisson impostors (no
clustering); SMF without redshift evolution; φ (population) and SMF⊗relation (host) not forced to be
one function; the selection horizon fitted on the table; timeouts / p0-window losses not in
P_det; completeness evaluated at the truth h (the advisor's FROZEN variant showed this is minor).
