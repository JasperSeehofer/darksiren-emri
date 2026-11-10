# BUDGET READOUT — Rank-1 information budget (M192 item 2)

**PROVISIONAL** — every number below is provisional until the author approves the derivation `docs/derivations/rank1_information_budget.md` (C-10; card `m192-r1-rank1` Q1). Pre-registration: `.planning/missions/M192/PREREGISTRATION_RANK1.md` (md5 `90dc4200f95f82c8a64b26f302f739c6`). Run: commit `14f7fd2365`, 2026-11-10T22:17:13.455896+00:00, seed base 192000, 8 realisations. Fisher, mock-internal.

σ_h is in units of h (truth h = 0.73). *cond* = conditional on the fiducial population, selection and relation; *marg* = marginalised over κ, s, a, σ_rel at the registered priors (κ 0.5, s 0.3, a 0.5 dex, σ_rel 0.2 dex). Information channels: *full* = host association + population; *population* = empty catalogue; *association* = the increment, (1/σ_full² − 1/σ_pop²)^(−1/2). Clustering: not computed.

## Kill criterion (H-MTC-popdom)

**R_GLADE = — ± nan** (seed61000, full information channel, all events, registered priors); seed62000 row: —. Line: 1.30.

## seed61000: all events

Calibrations: κ₀ = -1.35; P_det horizon μ = 1.360 + -0.107x + -1.307x² (x = ln(M_z/5e5)), σ_hor = 0.464; z_max = 1.44.

| scenario | information channel | σ_h cond | σ_h marg | R = marg/cond | σ_h/h cond | σ_h/h marg |
|---|---|---|---|---|---|---|
| GLADE+ today | full | 0.0500 | — | — ± nan | 6.85 % | — |
| GLADE+ today | population | 0.0500 | — | — | 6.86 % | — |
| GLADE+ today | association | 1.1615 | — | — | 159.11 % | — |
| LSST+Euclid photometry | full | 0.0493 | — | — ± nan | 6.76 % | — |
| LSST+Euclid photometry | population | 0.0500 | — | — | 6.86 % | — |
| LSST+Euclid photometry | association | 0.2946 | — | — | 40.36 % | — |
| pointed spectroscopy (N=20) | full | 0.0037 | 0.0037 | 1.006 ± 0.000 | 0.51 % | 0.51 % |
| pointed spectroscopy (N=20) | population | 0.0500 | — | — | 6.86 % | — |
| pointed spectroscopy (N=20) | association | 0.0037 | — | — | 0.51 % | — |

### By redshift bin (full information channel)

| scenario | z bin | N | σ_h cond | σ_h marg | R | population share of cond. information |
|---|---|---|---|---|---|---|
| GLADE+ today | z[0.0,0.25) | 174 | 0.0736 | 0.0803 | 1.092 | 0.996 |
| GLADE+ today | z[0.25,0.5) | 658 | 0.0543 | 0.0722 | 1.330 | 1.000 |
| GLADE+ today | z[0.5,0.75) | 614 | — | — | — | — |
| GLADE+ today | z[0.75,inf) | 144 | — | — | — | — |
| LSST+Euclid photometry | z[0.0,0.25) | 174 | 0.0731 | 0.0797 | 1.090 | 0.982 |
| LSST+Euclid photometry | z[0.25,0.5) | 658 | 0.0536 | 0.0706 | 1.318 | 0.975 |
| LSST+Euclid photometry | z[0.5,0.75) | 614 | — | — | — | — |
| LSST+Euclid photometry | z[0.75,inf) | 144 | — | — | — | — |
| pointed spectroscopy (N=20) | z[0.0,0.25) | 174 | 0.0037 | 0.0037 | 1.000 | 0.003 |
| pointed spectroscopy (N=20) | z[0.25,0.5) | 658 | 0.0543 | 0.0723 | 1.331 | 1.001 |
| pointed spectroscopy (N=20) | z[0.5,0.75) | 614 | — | — | — | — |
| pointed spectroscopy (N=20) | z[0.75,inf) | 144 | — | — | — | — |

### Catalogue census (mean over realisations; check, not calibrated)

| scenario | events with no candidate | median K (non-empty) | host catalogued share | hosts catalogued |
|---|---|---|---|---|
| GLADE+ today | 58 (3.6 %) | 613 | 0.01 % | 0.1 |
| LSST+Euclid photometry | 8 (0.5 %) | 31802 | 77.88 % | 1238.2 |
| pointed spectroscopy (N=20) | 0 (0.0 %) | 31440 | 77.70 % | 1235.5 |

Repo census for comparison (GLADE+, production 1D ball, 1,588 events): 606 empty, non-empty median K = 84; mock in-catalogue share 76/1,590 = 4.8 % (mass-independent Bernoulli in the generator).

### Post-hoc diagnostics (not registered): why the marginalised entries are undefined

*eig(F+P)* = eigenvalues of the registered observed-information matrix plus prior precision (order h, κ, s, a, σ_rel); a negative eigenvalue makes σ_h,marg undefined. *OPG* = the outer product of per-event scores, Σ_i s_i s_iᵀ, which is positive semi-definite by construction and equals the Fisher matrix in expectation only if the model is correctly specified; the gap between OPG and the observed information is a misspecification diagnostic (White 1982).

| scenario | min eig(F+P) | eigenvector of the min eigenvalue | σ_h cond (OPG) | σ_h marg (OPG) | R (OPG) | R (OPG), population only |
|---|---|---|---|---|---|---|
| GLADE+ today | -41.3 | h -0.42, kappa -0.90, s +0.08, a +0.00, sig_rel -0.00 | 0.0456 | 0.1243 | 2.725 | 2.726 |
| LSST+Euclid photometry | -39.4 | h -0.41, kappa -0.91, s +0.08, a +0.00, sig_rel +0.00 | 0.0459 | 0.1201 | 2.619 | 2.726 |
| pointed spectroscopy (N=20) | 4.06 | h +0.00, kappa +0.00, s +0.00, a -1.00, sig_rel +0.00 | 0.0037 | 0.0038 | 1.004 | 2.726 |

## Sensitivity and robustness rows (seed61000, all events, full information channel)

| row | scenario | σ_h cond | σ_h marg | R |
|---|---|---|---|---|
| priors_x0.5 | GLADE+ today | 0.0500 | — | — |
| priors_x0.5 | LSST+Euclid photometry | 0.0493 | — | — |
| priors_x0.5 | pointed spectroscopy (N=20) | 0.0037 | 0.0037 | 1.005 |
| priors_x2 | GLADE+ today | 0.0500 | — | — |
| priors_x2 | LSST+Euclid photometry | 0.0493 | — | — |
| priors_x2 | pointed spectroscopy (N=20) | 0.0037 | 0.0037 | 1.006 |
| kappa_prior_0.25 | GLADE+ today | 0.0500 | — | — |
| kappa_prior_0.25 | LSST+Euclid photometry | 0.0493 | — | — |
| kappa_prior_0.25 | pointed spectroscopy (N=20) | 0.0037 | 0.0037 | 1.005 |
| kappa_prior_1.0 | GLADE+ today | 0.0500 | — | — |
| kappa_prior_1.0 | LSST+Euclid photometry | 0.0493 | — | — |
| kappa_prior_1.0 | pointed spectroscopy (N=20) | 0.0037 | 0.0037 | 1.006 |
| pozzoli_N139 | GLADE+ today | 0.1654 | 1.0274 | 6.392 |
| pozzoli_N139 | LSST+Euclid photometry | 0.1647 | 0.6130 | 3.912 |
| pozzoli_N139 | pointed spectroscopy (N=20) | 0.0143 | 0.0137 | 1.004 |
| Om_free | GLADE+ today | 0.0488 | — | — |
| Om_free | LSST+Euclid photometry | 0.0482 | — | — |
| Om_free | pointed spectroscopy (N=20) | 0.0037 | 0.0038 | 1.015 |

## seed62000: all events

Calibrations: κ₀ = -1.15; P_det horizon μ = 1.254 + -0.422x + -0.890x² (x = ln(M_z/5e5)), σ_hor = 0.505; z_max = 1.32.

| scenario | information channel | σ_h cond | σ_h marg | R = marg/cond | σ_h/h cond | σ_h/h marg |
|---|---|---|---|---|---|---|
| GLADE+ today | full | 0.0552 | — | — ± nan | 7.57 % | — |
| GLADE+ today | population | 0.0552 | — | — | 7.57 % | — |
| GLADE+ today | association | 13.2133 | — | — | 1810.05 % | — |
| LSST+Euclid photometry | full | 0.0547 | — | — ± nan | 7.49 % | — |
| LSST+Euclid photometry | population | 0.0552 | — | — | 7.57 % | — |
| LSST+Euclid photometry | association | 0.3768 | — | — | 51.62 % | — |
| pointed spectroscopy (N=20) | full | 0.0031 | 0.0031 | 1.005 ± 0.000 | 0.42 % | 0.42 % |
| pointed spectroscopy (N=20) | population | 0.0552 | — | — | 7.57 % | — |
| pointed spectroscopy (N=20) | association | 0.0031 | — | — | 0.42 % | — |

### By redshift bin (full information channel)

| scenario | z bin | N | σ_h cond | σ_h marg | R | population share of cond. information |
|---|---|---|---|---|---|---|
| GLADE+ today | z[0.0,0.25) | 192 | 0.0662 | 0.0695 | 1.050 | 1.000 |
| GLADE+ today | z[0.25,0.5) | 586 | 0.0628 | 0.1244 | 1.980 | 1.000 |
| GLADE+ today | z[0.5,0.75) | 646 | — | — | — | — |
| GLADE+ today | z[0.75,inf) | 121 | — | — | — | — |
| LSST+Euclid photometry | z[0.0,0.25) | 192 | 0.0659 | 0.0691 | 1.049 | 0.990 |
| LSST+Euclid photometry | z[0.25,0.5) | 586 | 0.0621 | 0.1192 | 1.920 | 0.978 |
| LSST+Euclid photometry | z[0.5,0.75) | 646 | — | — | — | — |
| LSST+Euclid photometry | z[0.75,inf) | 121 | — | — | — | — |
| pointed spectroscopy (N=20) | z[0.0,0.25) | 192 | 0.0031 | 0.0031 | 1.000 | 0.002 |
| pointed spectroscopy (N=20) | z[0.25,0.5) | 586 | 0.0623 | 0.1203 | 1.932 | 0.983 |
| pointed spectroscopy (N=20) | z[0.5,0.75) | 646 | — | — | — | — |
| pointed spectroscopy (N=20) | z[0.75,inf) | 121 | — | — | — | — |

### Catalogue census (mean over realisations; check, not calibrated)

| scenario | events with no candidate | median K (non-empty) | host catalogued share | hosts catalogued |
|---|---|---|---|---|
| GLADE+ today | 54 (3.5 %) | 593 | 0.02 % | 0.2 |
| LSST+Euclid photometry | 2 (0.2 %) | 31524 | 77.50 % | 1197.4 |
| pointed spectroscopy (N=20) | 0 (0.0 %) | 31452 | 77.86 % | 1202.9 |

Repo census for comparison (GLADE+, production 1D ball, 1,588 events): 606 empty, non-empty median K = 84; mock in-catalogue share 76/1,590 = 4.8 % (mass-independent Bernoulli in the generator).

### Post-hoc diagnostics (not registered): why the marginalised entries are undefined

*eig(F+P)* = eigenvalues of the registered observed-information matrix plus prior precision (order h, κ, s, a, σ_rel); a negative eigenvalue makes σ_h,marg undefined. *OPG* = the outer product of per-event scores, Σ_i s_i s_iᵀ, which is positive semi-definite by construction and equals the Fisher matrix in expectation only if the model is correctly specified; the gap between OPG and the observed information is a misspecification diagnostic (White 1982).

| scenario | min eig(F+P) | eigenvector of the min eigenvalue | σ_h cond (OPG) | σ_h marg (OPG) | R (OPG) | R (OPG), population only |
|---|---|---|---|---|---|---|
| GLADE+ today | -73.5 | h +0.49, kappa +0.87, s -0.09, a +0.00, sig_rel -0.00 | 0.0468 | 0.1295 | 2.766 | 2.766 |
| LSST+Euclid photometry | -71.8 | h +0.49, kappa +0.87, s -0.09, a -0.00, sig_rel -0.00 | 0.0467 | 0.1243 | 2.662 | 2.766 |
| pointed spectroscopy (N=20) | 4.16 | h +0.00, kappa +0.00, s -0.00, a -1.00, sig_rel +0.01 | 0.0028 | 0.0028 | 1.002 | 2.766 |

## Verdict and deviations from the registration

**Registered verdict: INCONCLUSIVE** (PREREGISTRATION_RANK1.md §5: 'INCONCLUSIVE ... if F + P is not positive definite'). In the GLADE+ and LSST scenarios the observed-information matrix plus prior precision has one negative eigenvalue, along the (h, κ) direction, on both tables. σ_h,marg and so R_GLADE are undefined at the registered setting. The post-hoc OPG diagnostic, which is positive semi-definite by construction, gives R ≈ 2.7 on both tables: well above 1.30, and population-only ≈ full. It is not registered and does not decide the kill criterion. The Pozzoli-like robustness row (N = 139) is positive definite with R_GLADE = 6.4.

Numerical deviations from §1 (implementation, not model): (i) a per-event local z grid (161 nodes over ±10 σ_z of the GW likelihood + 1 % of z) replaces the registered global 2,000-node grid. A first run with a ±(9σ_z + 0.003(1+z)) window under-resolved the events with σ_dL/d_L ≈ 1e-3 (observed curvature −2.3e4 for one event) and was discarded before any number was read against the criterion. (ii) The census K counts zhat inside the z range of d_L ± 3σ over h ∈ [0.60, 0.86], the production ball's analogue. (iii) The rows for single redshift bins carry the whole-sample α curvature per event, as registered. Observed information per bin is negative above z = 0.5, so those rows are blank: the model's p(z) and P_det do not reproduce the table's high-z tail (misspecification, see the diagnostics).
