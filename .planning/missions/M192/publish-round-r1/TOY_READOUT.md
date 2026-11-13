# TOY READOUT — line-of-sight toy with a mass dimension (M192 item 4)

Registered design: `.planning/missions/M192/PREREGISTRATION_RANK1.md` §6 (md5 `90dc4200…`). Code: `analysis/rank1/toy/q4_mass_toy.py`; run metadata `run_metadata.json` (seeds 4192001–4192072, one per universe). 72 universes × 300 events, truth h = 0.73, h grid 0.43–1.03 (121 nodes). In-catalogue share of detected hosts: 0.020.

Pulls are (posterior mean − **truth 0.73**)/σ_h (C-21). s.e. of the median pull = 1.2533·sd/√72. *scaled* = median pull × √(1588/300) (a fixed bias with σ ∝ N^−½). *MAD* = 1.4826 × median absolute deviation of the pulls (robust spread; sd is inflated by a few spuriously narrow posteriors).

| variant | ⟨σ_h⟩ | ⟨mean⟩ − truth | median pull ± s.e. | MAD pull | cov68 / cov95 | ⟨PIT⟩ | median pull scaled to N = 1,588 |
|---|---|---|---|---|---|---|---|
| CORRECT | 0.057 | +0.011 | +0.12 ± 0.15 | 0.80 | 0.69 / 0.93 | 0.53 | +0.28 |
| FINITE-α (4.8×) | 0.043 | -0.129 | -2.80 ± 3.22 | 5.64 | 0.21 / 0.36 | 0.78 | -6.45 |
| FINITE-α (9.6×) | 0.043 | -0.129 | -2.65 ± 3.22 | 4.25 | 0.25 / 0.46 | 0.78 | -6.11 |
| DL-CUT-NOT-IN-α | 0.047 | +0.020 | +0.23 ± 0.16 | 0.86 | 0.67 / 0.94 | 0.44 | +0.53 |
| MASS-SEL-NOT-IN-α | 0.056 | -0.055 | -0.39 ± 0.27 | 0.73 | 0.67 / 0.90 | 0.67 | -0.89 |
| MIX (×2.3) | 0.049 | +0.029 | +0.39 ± 0.20 | 1.03 | 0.56 / 0.90 | 0.47 | +0.90 |

**Post-hoc diagnostics — not registered.** Added after a development run on non-registered seeds (900–907) showed the registered FINITE variants biased far beyond pool noise. They decompose FINITE-α into its two parts:

| diagnostic | ⟨σ_h⟩ | ⟨mean⟩ − truth | median pull ± s.e. | MAD pull | cov68 / cov95 | ⟨PIT⟩ | scaled |
|---|---|---|---|---|---|---|---|
| INF-POOL: α = exact population expectation W·f(h) (infinite pool) | 0.061 | -0.136 | -1.30 ± 0.61 | 1.74 | 0.26 / 0.57 | 0.84 | -3.00 |
| FINITE-DARK (4.8×): catalogue part of α exact, dark part from the pool | 0.052 | -0.032 | -0.31 ± 1.63 | 2.09 | 0.32 / 0.57 | 0.58 | -0.71 |
| FINITE-DARK (9.6×) | 0.061 | -0.019 | -0.01 ± 0.45 | 1.25 | 0.47 / 0.74 | 0.54 | -0.02 |

## Read against the registration (§6)

| registered prediction | outcome | verdict |
|---|---|---|
| CORRECT closes, \|median pull\| < 0.25 | +0.12 ± 0.15; cov68 0.69, cov95 0.93 | **met** |
| MIX positive, ≈ +0.5 | +0.39 ± 0.20 | **met** (wrong sign for the repo) |
| FINITE-α: no systematic sign, extra scatter, halved at 9.6× | median -2.80 (4.8×), -2.65 (9.6×); ⟨mean⟩ − truth −0.129 in both | **not met** — a large systematic negative offset that does not shrink with the pool |
| MASS-SEL-NOT-IN-α negative, ≈ −1 (≈ −2.3 scaled) | -0.39 ± 0.27 (-0.89 scaled) | sign met, size **not met** (about 40 %) |
| DL-CUT-NOT-IN-α negative, \|pull\| ≈ 0.5 | +0.23 ± 0.16 | **not met** (positive, within 1.5 s.e. of 0) |
| no single variant reaches −3.5σ at the repo's scale | registered FINITE-α reaches −6 (scaled), with very broad spread | **not met** — FINITE-α alone is sufficient in sign and size |

## What it says about the offset

1. **The selection explanation that produces a negative offset of the repo's size is the normaliser, not the pool size.** The registered FINITE-α variants take the whole α(h) from an injection pool drawn from the population. They shift h by −0.129 (≈ −2.7 median pull, −6 scaled) at both pool sizes. The diagnostics separate the parts. An infinite pool (INF-POOL, the exact population expectation) carries the same −0.136. With the catalogue part of α exact and only the dark part from the pool (FINITE-DARK), the median offset is −0.31 at 4.8× and −0.01 at 9.6×. The pool-noise pathology roughly halves with the doubled pool (cov68 0.32 → 0.47, MAD 2.1 → 1.25).
   **Mechanism.** In a Gray-type likelihood the catalogue term sums the realised catalogue at the true h, while the completion term uses the completeness at trial h. α must be built with the same split. A population-level α cannot see the realised catalogue and is mismatched. Its h-slope error is a few % per unit h, and multiplied by N σ_h² it gives a shift of order −0.1. This is the advisor's ranking #1 (leg-inconsistent normalisation) made concrete. Whether the repo's α has this structure is a round-2 question, not established here.
2. The mass-dependent loss absent from α (timeouts/p0 window) gives the repo's sign at about −0.9σ scaled to 1,588 events. It is a contributor, not sufficient alone.
3. The σ_dL/d_L cut and the ×2.3 mixture weight give positive pulls: they cannot explain a negative offset.
4. With deterministic (oracle) selection a finite pool makes α(h) a step function in h. At 4.8× this produces spuriously narrow posteriors (sd of pull ≫ MAD). This is a toy artefact that a smooth P_det would soften.

**Limits.** 1-D patches, Poisson galaxies, exact spec-z, oracle selection, toy horizon; 300 events per universe. A toy localises mechanisms; it does not reproduce the repo's numbers.
