"""Render results/rank1_<date>/budget.json as BUDGET_READOUT.md (GLADE+ today first in every table).

usage: python -m analysis.rank1.readout results/rank1_20261110
"""

import json
import math
import sys
from pathlib import Path
from typing import Any

KEYS = ["GLADE+ today", "LSST+Euclid photometry", "pointed spectroscopy (N=20)"]
H_TRUE = 0.73


def f(x: Any, nd: int = 4) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{x:.{nd}f}"


def pct(x: Any) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "—"
    return f"{100 * x / H_TRUE:.2f} %"


def channel_table(reg: dict[str, Any], mask: str) -> list[str]:
    out = [
        "| scenario | information channel | σ_h cond | σ_h marg | R = marg/cond | σ_h/h cond | σ_h/h marg |",
        "|---|---|---|---|---|---|---|",
    ]
    for k in KEYS:
        r = reg[k][mask]
        for ch in ("FULL", "POPULATION", "ASSOCIATION"):
            c = r[ch]
            R = c.get("R")
            Rs = f(R, 3) + (f" ± {c['R_se']:.3f}" if ch == "FULL" and "R_se" in c else "")
            out.append(
                f"| {k} | {ch.lower()} | {f(c['cond'])} | {f(c['marg'])} | {Rs if R is not None else '—'} "
                f"| {pct(c['cond'])} | {pct(c['marg'])} |"
            )
    return out


def main(argv: list[str]) -> int:
    d = Path(argv[1])
    J = json.loads((d / "budget.json").read_text())
    L: list[str] = []
    L.append("# BUDGET READOUT — Rank-1 information budget (M192 item 2)")
    L.append("")
    L.append(
        "**PROVISIONAL** — every number below is provisional until the author approves the derivation "
        "`docs/derivations/rank1_information_budget.md` (C-10; card `m192-r1-rank1` Q1). Pre-registration: "
        "`.planning/missions/M192/PREREGISTRATION_RANK1.md` (md5 `90dc4200f95f82c8a64b26f302f739c6`). "
        f"Run: commit `{J['meta']['git_commit'][:10]}`, {J['meta']['timestamp']}, seed base "
        f"{J['meta']['seed_base']}, {J['meta']['n_real']} realisations. Fisher, mock-internal."
    )
    L.append("")
    L.append(
        "σ_h is in units of h (truth h = 0.73). *cond* = conditional on the fiducial population, selection and "
        "relation; *marg* = marginalised over κ, s, a, σ_rel at the registered priors (κ 0.5, s 0.3, a 0.5 dex, "
        "σ_rel 0.2 dex). Information channels: *full* = host association + population; *population* = empty "
        "catalogue; *association* = the increment, (1/σ_full² − 1/σ_pop²)^(−1/2). Clustering: not computed."
    )
    s61 = J["seed61000"]
    R = s61["registered"]["GLADE+ today"]["all"]["FULL"]
    L.append("")
    L.append("## Kill criterion (H-MTC-popdom)")
    L.append("")
    s62R = (
        J.get("seed62000", {})
        .get("registered", {})
        .get("GLADE+ today", {})
        .get("all", {})
        .get("FULL", {})
        .get("R")
    )
    L.append(
        f"**R_GLADE = {f(R['R'], 3)} ± {R['R_se']:.3f}** (seed61000, full information channel, all events, "
        f"registered priors); seed62000 row: {f(s62R, 3)}. Line: 1.30."
    )
    for tname in ("seed61000", "seed62000"):
        if tname not in J:
            continue
        T = J[tname]
        L.append("")
        L.append(f"## {tname}: all events")
        L.append("")
        L.append(
            f"Calibrations: κ₀ = {T['kappa0']:.2f}; P_det horizon μ = {T['selection']['c0']:.3f} + "
            f"{T['selection']['c1']:.3f}x + {T['selection']['c2']:.3f}x² (x = ln(M_z/5e5)), σ_hor = "
            f"{T['selection']['sig']:.3f}; z_max = {T['zmax']:.2f}."
        )
        L.append("")
        L.extend(channel_table(T["registered"], "all"))
        L.append("")
        L.append("### By redshift bin (full information channel)")
        L.append("")
        L.append(
            "| scenario | z bin | N | σ_h cond | σ_h marg | R | population share of cond. information |"
        )
        L.append("|---|---|---|---|---|---|---|")
        for k in KEYS:
            for m, r in T["registered"][k].items():
                if m == "all" or m.endswith("__diag"):
                    continue
                c = r["FULL"]
                L.append(
                    f"| {k} | {m} | {r['n_events']} | {f(c['cond'])} | {f(c['marg'])} | {f(c['R'], 3)} | "
                    f"{f(r['pop_share_cond_info'], 3)} |"
                )
        L.append("")
        L.append("### Catalogue census (mean over realisations; check, not calibrated)")
        L.append("")
        L.append(
            "| scenario | events with no candidate | median K (non-empty) | host catalogued share | hosts catalogued |"
        )
        L.append("|---|---|---|---|---|")
        for k in KEYS:
            c = T["census_mean"][k]
            L.append(
                f"| {k} | {c['n_empty']:.0f} ({100 * c['frac_empty']:.1f} %) | {c['median_K_nonempty']:.0f} | "
                f"{100 * c['host_catalogued_share']:.2f} % | {c.get('host_catalogued_n', float('nan')):.1f} |"
            )
        L.append("")
        L.append(
            "Repo census for comparison (GLADE+, production 1D ball, 1,588 events): 606 empty, non-empty median "
            "K = 84; mock in-catalogue share 76/1,590 = 4.8 % (mass-independent Bernoulli in the generator)."
        )
        L.append("")
        L.append(
            "### Post-hoc diagnostics (not registered): why the marginalised entries are undefined"
        )
        L.append("")
        L.append(
            "*eig(F+P)* = eigenvalues of the registered observed-information matrix plus prior precision "
            "(order h, κ, s, a, σ_rel); a negative eigenvalue makes σ_h,marg undefined. *OPG* = the "
            "outer product of per-event scores, Σ_i s_i s_iᵀ, which is positive semi-definite by construction "
            "and equals the Fisher matrix in expectation only if the model is correctly specified; the gap "
            "between OPG and the observed information is a misspecification diagnostic (White 1982)."
        )
        L.append("")
        L.append(
            "| scenario | min eig(F+P) | eigenvector of the min eigenvalue | σ_h cond (OPG) | σ_h marg (OPG) | R (OPG) | R (OPG), population only |"
        )
        L.append("|---|---|---|---|---|---|---|")
        for k in KEYS:
            dg = T["registered"][k]["all__diag"]
            vec = ", ".join(f"{a} {b:+.2f}" for a, b in dg["eigvec_min"].items())
            o = dg["FULL_OPG"]
            po = dg["POPULATION_OPG"]
            L.append(
                f"| {k} | {dg['eig_F_plus_P'][0]:.3g} | {vec} | {f(o['cond'])} | {f(o['marg'])} | {f(o['R'], 3)} | {f(po['R'], 3)} |"
            )
        if "sensitivity" in T:
            L.append("")
            L.append(
                "## Sensitivity and robustness rows (seed61000, all events, full information channel)"
            )
            L.append("")
            L.append("| row | scenario | σ_h cond | σ_h marg | R |")
            L.append("|---|---|---|---|---|")
            for row, block in T["sensitivity"].items():
                for k in KEYS:
                    c = (
                        block[k]["N139"]["FULL"]
                        if row != "pozzoli_N139" and "N139" in block[k]
                        else None
                    )
                    if row == "pozzoli_N139":
                        c = block[k]["FULL"]
                    elif c is None:
                        c = block[k]["all"]["FULL"]
                    L.append(
                        f"| {row} | {k} | {f(c['cond'])} | {f(c['marg'])} | {f(c.get('R'), 3)} |"
                    )
    L.append("")
    L.append("## Verdict and deviations from the registration")
    L.append("")
    L.append(
        "**Registered verdict: INCONCLUSIVE** (PREREGISTRATION_RANK1.md §5: 'INCONCLUSIVE ... if F + P is not "
        "positive definite'). In the GLADE+ and LSST scenarios the observed-information matrix plus prior "
        "precision has one negative eigenvalue, along the (h, κ) direction, on both tables. σ_h,marg and so "
        "R_GLADE are undefined at the registered setting. The post-hoc OPG diagnostic, which is positive "
        "semi-definite by construction, gives R ≈ 2.7 on both tables: well above 1.30, and population-only "
        "≈ full. It is not registered and does not decide the kill criterion. The Pozzoli-like robustness "
        "row (N = 139) is positive definite with R_GLADE = 6.4."
    )
    L.append("")
    L.append(
        "Numerical deviations from §1 (implementation, not model): (i) a per-event local z grid (161 nodes "
        "over ±10 σ_z of the GW likelihood + 1 % of z) replaces the registered global 2,000-node grid. A first "
        "run with a ±(9σ_z + 0.003(1+z)) window under-resolved the events with σ_dL/d_L ≈ 1e-3 "
        "(observed curvature −2.3e4 for one event) and was discarded before any number was read against the "
        "criterion. (ii) The census K counts zhat inside the z range of d_L ± 3σ over h ∈ [0.60, 0.86], "
        "the production ball's analogue. (iii) The rows for single redshift bins carry the whole-sample α "
        "curvature per event, as registered. Observed information per bin is negative above z = 0.5, so "
        "those rows are blank: the model's p(z) and P_det do not reproduce the table's high-z tail "
        "(misspecification, see the diagnostics)."
    )
    (d / "BUDGET_READOUT.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
