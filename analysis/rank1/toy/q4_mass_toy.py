"""Line-of-sight toy with a mass dimension (mission M192, scope item 4).

Extends the outside advisor's q4_line_toy.py (verbatim copy next to this file; garden commit
55a3be19, md5 d3e8c73c) with:
  * galaxy stellar masses from the Baldry+12 SMF (luminosity proportional to M*; flux-limited
    catalogue x = M*/10^10.66 > FLIM d_L(z, h_true)^2, as the base toy's L/L*);
  * M_BH ~ lognormal(a + beta (log M* - 11), sigma_rel), a = 7.45, beta = 1.05, sigma = 0.55;
  * hosts drawn with weight w(M) = M^-0.13 on [1e4, 1e7] M_sun (rate per BH), among galaxies whose
    event is detectable: d_L(z, h_true) < D_hor(M_z), D_hor lognormal in M_z (as the q3 toy);
  * estimator: catalogue galaxy g weighs p(log Mhat_z - log(1+z_g) | M*_g) w(.) and the dark term
    uses Phi_dark = int n(M*) [1 - c(M*, z, h)] p(M|M*) w(M) dlog M*; selection in numerator and
    alpha (theta-deterministic);
  * six variants (PREREGISTRATION_RANK1.md section 6): CORRECT, FINITE-alpha (4.8x, 9.6x, pool
    drawn inside the toy), DL-CUT-NOT-IN-alpha, MASS-SEL-NOT-IN-alpha, MIX (x2.3).

usage: python -m analysis.rank1.toy.q4_mass_toy --seed0 4192001 --n_universes 72 --out DIR
Every universe u uses seed seed0 + u (recorded); run_metadata.json is written in DIR (C-16).
"""

import argparse
import datetime as _dt
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt

Arr = npt.NDArray[np.float64]

OM, C, HT, SIG, NEV, P, NBOX = 0.2726, 299792.458, 0.73, 0.04, 300, 600, 1000.0
FLIM = 6e-7
A_REL, B_REL, S_REL = 7.45, 1.05, 0.55
VARIANTS = ["CORRECT", "FINITE4.8", "FINITE9.6", "DLCUT", "MASSSEL", "MIX2.3"]
# Post-hoc DIAGNOSTICS (not registered; added after a dev run on non-registered seeds showed that
# a population-drawn alpha cannot see the realised catalogue): INFPOOL = exact population
# expectation W*f(h) for the whole alpha; FINITE*-DARK = catalogue part of alpha exact (from the
# realised catalogue), only the dark part from the finite pool.
DIAGNOSTICS = ["INFPOOL", "FINITE4.8-DARK", "FINITE9.6-DARK"]
MIX_W = 2.3

# grids
ZG = np.linspace(0.0, 0.6, 3001)
_EZ = np.sqrt(OM * (1 + ZG) ** 3 + 1 - OM)
_DC = np.concatenate([[0.0], np.cumsum(0.5 * (1 / _EZ[1:] + 1 / _EZ[:-1]) * np.diff(ZG))]) * C / 100
ZG, _EZ, _DC = ZG[1:], _EZ[1:], _DC[1:]
_DVDZ = _DC**2 / _EZ
NB_FINE = NBOX * _DVDZ / np.trapezoid(_DVDZ, ZG)
ZE = np.linspace(ZG[0], 0.6, 500)  # estimator z grid
NB_E = np.interp(ZE, ZG, NB_FINE)
LMS = np.linspace(8.0, 12.0, 81)  # log10 M*
LM = np.linspace(4.0, 7.0, 61)  # log10 M_BH
DLM = LM[1] - LM[0]
HGRID = np.linspace(0.43, 1.03, 121)


def dLz(z: Arr, h: float) -> Arr:
    return np.asarray((1 + z) * np.interp(z, ZG, _DC) / h)


DH0 = float(dLz(np.array([0.45]), HT)[0])


def d_hor(lmz: Arr) -> Arr:
    """Detection horizon (Mpc) vs log10 detector-frame mass: peak at 10^5.8, width 0.6 dex."""
    return np.asarray(DH0 * np.exp(-0.5 * ((lmz - 5.8) / 0.6) ** 2))


def smf_w(lms: Arr) -> Arr:
    x = 10 ** (lms - 10.66)
    return np.asarray(np.exp(-x) * (3.96e-3 * x ** (-0.35 + 1) + 0.79e-3 * x ** (-1.47 + 1)))


NFRAC = smf_w(LMS) / smf_w(LMS).sum()
XMS = 10 ** (LMS - 10.66)


def w_rate(lm: Arr) -> Arr:
    return np.where((lm >= 4.0) & (lm <= 7.0), 10 ** (-0.13 * (lm - 5.0)), 0.0)


def p_ln(lm: Arr, lms: Arr) -> Arr:
    mu = A_REL + B_REL * (lms - 11.0)
    return np.asarray(np.exp(-0.5 * ((lm - mu) / S_REL) ** 2) / (np.sqrt(2 * np.pi) * S_REL))


# ----------------------------------------------------------------------------- universe
def sky(rng: np.random.Generator) -> dict[str, Any]:
    n = rng.poisson(NBOX, P)
    patch = np.repeat(np.arange(P), n)
    tot = int(n.sum())
    z = np.abs(rng.choice(ZG, tot, p=NB_FINE / NB_FINE.sum()) + rng.uniform(-1e-4, 1e-4, tot))
    k = rng.choice(len(LMS), tot, p=NFRAC)
    lms = LMS[k]
    lm = A_REL + B_REL * (lms - 11.0) + S_REL * rng.standard_normal(tot)
    incat = XMS[k] > FLIM * dLz(z, HT) ** 2
    q = np.exp(0.5 * rng.standard_normal(tot))  # DL-CUT degeneracy factor
    return {"z": z, "patch": patch, "k": k, "lm": lm, "incat": incat, "q": q}


def draw_events(s: dict[str, Any], rng: np.random.Generator, variant: str) -> npt.NDArray[np.int64]:
    z, lm = s["z"], s["lm"]
    lmz = lm + np.log10(1 + z)
    dl = dLz(z, HT)
    ok = dl < d_hor(lmz)
    if variant == "DLCUT":
        r = 0.05 * dl / d_hor(lmz) * s["q"]
        ok &= r < 0.10
    if variant == "MASSSEL":
        ok &= (lm >= np.log10(2e5)) & (lm <= np.log10(3.3e6))
    w = w_rate(lm) * ok
    idx = np.where(w > 0)[0]
    return np.asarray(rng.choice(idx, NEV, replace=False, p=w[idx] / w[idx].sum()))


# ----------------------------------------------------------------------------- likelihood
def numerators(
    s: dict[str, Any], hosts: npt.NDArray[np.int64], rng: np.random.Generator
) -> tuple[Arr, Arr]:
    """Catalogue and dark numerator terms on HGRID: (cat (NH, NEV), dark (NH, NEV))."""
    z, patch, lm = s["z"], s["patch"], s["lm"]
    zt = z[hosts]
    lmz = lm[hosts] + np.log10(1 + zt)
    dhat = dLz(zt, HT) * (1 + SIG * rng.standard_normal(NEV))
    dh = d_hor(lmz)
    # catalogue galaxies by patch
    cat = np.where(s["incat"])[0]
    pc = patch[cat]
    order = np.argsort(pc)
    cat = cat[order]
    pc = pc[order]
    hp = patch[hosts]
    st = np.searchsorted(pc, hp)
    en = np.searchsorted(pc, hp, side="right")
    # dark: G[i, z, k] cumulative over M* of nfrac * p_ln * w at the implied source mass
    lsrc = lmz[:, None] - np.log10(1 + ZE)[None, :]  # (NEV, NZ)
    integ = (
        NFRAC[None, None, :] * p_ln(lsrc[:, :, None], LMS[None, None, :]) * w_rate(lsrc)[:, :, None]
    )
    G = np.concatenate([np.zeros((NEV, len(ZE), 1)), np.cumsum(integ, axis=2)], axis=2)
    catT = np.zeros((len(HGRID), NEV))
    darkT = np.zeros((len(HGRID), NEV))
    cat_mass: list[Arr] = []
    cat_z: list[Arr] = []
    for i in range(NEV):
        g = cat[st[i] : en[i]]
        zg = z[g]
        ls = lmz[i] - np.log10(1 + zg)
        cat_mass.append(p_ln(ls, LMS[s["k"][g]]) * w_rate(ls))
        cat_z.append(zg)
    for j, h in enumerate(HGRID):
        dl = dLz(ZE, h)
        xlim = FLIM * dLz(ZE, h) ** 2
        kidx = np.searchsorted(XMS, xlim)  # number of M* cells with x < xlim (uncatalogued)
        Gz = np.take_along_axis(G, np.broadcast_to(kidx[None, :, None], (NEV, len(ZE), 1)), axis=2)[
            :, :, 0
        ]
        sel = dl[None, :] < dh[:, None]
        gw = np.exp(-0.5 * ((dhat[:, None] - dl[None, :]) / (SIG * dl[None, :])) ** 2) / dl[None, :]
        darkT[j] = np.trapezoid(NB_E[None, :] * Gz * sel * gw, ZE, axis=1)
        for i in range(NEV):
            if len(cat_z[i]) == 0:
                continue
            dlc = dLz(cat_z[i], h)
            gi = np.exp(-0.5 * ((dhat[i] - dlc) / (SIG * dlc)) ** 2) / dlc
            catT[j, i] = np.sum(cat_mass[i] * gi * (dlc < dh[i]))
    return catT, darkT


def alpha_exact(s: dict[str, Any]) -> tuple[Arr, Arr]:
    """(cat part per patch, dark part) of alpha on HGRID."""
    cat = np.where(s["incat"])[0]
    zc = s["z"][cat]
    A = p_ln(LM[None, :], LMS[s["k"][cat]][:, None]) * w_rate(LM)[None, :] * DLM  # (Ncat, NM)
    lmz_c = LM[None, :] + np.log10(1 + zc)[:, None]
    dhc = d_hor(lmz_c)
    B = NFRAC[:, None] * p_ln(LM[None, :], LMS[:, None]) * w_rate(LM)[None, :] * DLM  # (NMS, NM)
    lmz_d = LM[None, :] + np.log10(1 + ZE)[:, None]
    dhd = d_hor(lmz_d)  # (NZ, NM)
    Bc = np.concatenate([np.zeros((1, len(LM))), np.cumsum(B, axis=0)], axis=0)
    ac = np.zeros(len(HGRID))
    ad = np.zeros(len(HGRID))
    for j, h in enumerate(HGRID):
        ac[j] = np.sum(A * (dLz(zc, h)[:, None] < dhc)) / P
        dl = dLz(ZE, h)
        kidx = np.searchsorted(XMS, FLIM * dl**2)
        Bz = Bc[kidx]  # (NZ, NM)
        ad[j] = np.trapezoid(NB_E * np.sum(Bz * (dl[:, None] < dhd), axis=1), ZE)
    return ac, ad


def pop_total_and_fdet() -> tuple[float, Arr]:
    """Total population weight W and the exact detected fraction on HGRID."""
    B = NFRAC[:, None] * p_ln(LM[None, :], LMS[:, None]) * w_rate(LM)[None, :] * DLM
    bm = B.sum(axis=0)  # (NM,)
    lmz = LM[None, :] + np.log10(1 + ZE)[:, None]
    dh = d_hor(lmz)
    W = float(np.trapezoid(NB_E * bm.sum(), ZE))
    f = np.array(
        [
            np.trapezoid(NB_E * np.sum(bm[None, :] * (dLz(ZE, h)[:, None] < dh), axis=1), ZE) / W
            for h in HGRID
        ]
    )
    return W, f


def alpha_finite(rng: np.random.Generator, ratio: float, W: float, fdet: Arr) -> tuple[Arr, Arr]:
    """Finite injection pool drawn from the population: (total alpha, dark part of alpha)."""
    jt = int(np.argmin(np.abs(HGRID - HT)))
    n_inj = int(round(ratio * NEV / fdet[jt]))
    zi = rng.choice(ZG, n_inj, p=NB_FINE / NB_FINE.sum())
    # source-mass draw from sum_M* nfrac p_ln(.|M*) w(.) via rejection on w
    k = rng.choice(len(LMS), 6 * n_inj, p=NFRAC)
    lm = A_REL + B_REL * (LMS[k] - 11.0) + S_REL * rng.standard_normal(6 * n_inj)
    acc = rng.random(6 * n_inj) < w_rate(lm) / 10 ** (-0.13 * (4.0 - 5.0))
    lm = lm[acc][:n_inj]
    kk = k[acc][:n_inj]
    if len(lm) < n_inj:
        raise RuntimeError("injection rejection sampler ran short")
    dh = d_hor(lm + np.log10(1 + zi))
    tot = np.zeros(len(HGRID))
    dark = np.zeros(len(HGRID))
    for j, h in enumerate(HGRID):
        dl = dLz(zi, h)
        det = dl < dh
        tot[j] = W * np.mean(det)
        dark[j] = W * np.mean(det & (XMS[kk] < FLIM * dl**2))
    return tot, dark


def summ(LL: Arr) -> tuple[float, float, float, float, float]:
    p = np.exp(LL - LL.max())
    p /= p.sum()
    m = float((p * HGRID).sum())
    cdf = np.cumsum(p)
    sd = float(np.sqrt((p * (HGRID - m) ** 2).sum()))
    return (
        m,
        sd,
        float(np.interp(HT, HGRID, cdf)),
        float(p[0] + p[-1]),
        float(np.interp(0.5, cdf, HGRID)),
    )


def universe(seed: int, W: float, fdet: Arr) -> dict[str, Any]:
    # independent streams per component, so that variants do not shift each other's draws
    r_sky, r_ev, r_num, r_pool, r_dl, r_ms = (
        np.random.default_rng(c) for c in np.random.SeedSequence(seed).spawn(6)
    )
    s = sky(r_sky)
    ac, ad = alpha_exact(s)
    out: dict[str, Any] = {"seed": seed}
    base = draw_events(s, r_ev, "CORRECT")
    out["share"] = float(s["incat"][base].mean())
    catT, darkT = numerators(s, base, r_num)

    def ll(cw: float, alpha: Arr, cT: Arr, dT: Arr) -> Arr:
        return np.asarray(np.log(cw * cT + dT + 1e-300).sum(axis=1) - NEV * np.log(alpha))

    out["CORRECT"] = summ(ll(1.0, ac + ad, catT, darkT))
    out["MIX2.3"] = summ(ll(MIX_W, MIX_W * ac + ad, catT, darkT))
    for ratio in (4.8, 9.6):
        tot, dark = alpha_finite(r_pool, ratio, W, fdet)
        out[f"FINITE{ratio}"] = summ(ll(1.0, tot, catT, darkT))
        out[f"FINITE{ratio}-DARK"] = summ(ll(1.0, ac + dark, catT, darkT))
    out["INFPOOL"] = summ(ll(1.0, W * fdet, catT, darkT))
    for v, rv in (("DLCUT", r_dl), ("MASSSEL", r_ms)):
        hv = draw_events(s, rv, v)
        cT, dT = numerators(s, hv, rv)
        out[v] = summ(ll(1.0, ac + ad, cT, dT))
    return out


def readout(res: list[dict[str, Any]]) -> dict[str, Any]:
    table: dict[str, Any] = {}
    for v in VARIANTS + DIAGNOSTICS:
        r = np.array([u[v] for u in res])
        pull = (r[:, 0] - HT) / r[:, 1]
        n = len(pull)
        pit = r[:, 2]
        table[v] = {
            "n_universes": n,
            "mean_sigma_h": float(r[:, 1].mean()),
            "mean_minus_truth": float(r[:, 0].mean() - HT),
            "median_pull": float(np.median(pull)),
            "median_pull_se": float(1.2533 * pull.std(ddof=1) / np.sqrt(n)),
            "mean_pull": float(pull.mean()),
            "rms_pull": float(np.sqrt(np.mean(pull**2))),
            "sd_pull": float(pull.std(ddof=1)),
            "mad_pull": float(1.4826 * np.median(np.abs(pull - np.median(pull)))),
            "cov68": float(np.mean((pit > 0.16) & (pit < 0.84))),
            "cov95": float(np.mean((pit > 0.025) & (pit < 0.975))),
            "mean_PIT": float(pit.mean()),
            "median_pull_scaled_N1588": float(np.median(pull) * np.sqrt(1588 / NEV)),
            "mean_edge_mass": float(r[:, 3].mean()),
        }
    table["in_catalogue_share"] = float(np.mean([u["share"] for u in res]))
    return table


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed0", type=int, required=True)
    ap.add_argument("--n_universes", type=int, default=72)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    meta = {
        "git_commit": commit,
        "timestamp": _dt.datetime.now(_dt.UTC).isoformat(),
        "argv": sys.argv,
        "seeds": [args.seed0 + u for u in range(args.n_universes)],
        "prereg": ".planning/missions/M192/PREREGISTRATION_RANK1.md",
    }
    (args.out / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    W, fdet = pop_total_and_fdet()
    res = []
    t0 = time.time()
    for u in range(args.n_universes):
        res.append(universe(args.seed0 + u, W, fdet))
        print(
            f"universe {u + 1}/{args.n_universes} seed {args.seed0 + u} {time.time() - t0:.0f}s",
            flush=True,
        )
        (args.out / "toy_universes.json").write_text(json.dumps(res))
    tab = readout(res)
    (args.out / "toy_readout.json").write_text(json.dumps(tab, indent=2))
    print(json.dumps(tab, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
