"""Rank-1 information budget (mission M192, scope item 2).

Fisher budget of the H0 information in the frozen EMRI event table, conditional on and
marginalised over the population (rate evolution kappa, mass-function slope s), and the
black-hole--galaxy relation (normalisation a, scatter sigma_rel), per information channel
(full = host association + population; population only = empty catalogue; association
increment), redshift bin and catalogue scenario (GLADE+ today first).

Model and derivation: docs/derivations/rank1_information_budget.md (C-10: presented /
PENDING-AUTHOR; every number is provisional until the author approves it).
Pre-registration: .planning/missions/M192/PREREGISTRATION_RANK1.md (C-14).

Read-only imports from darksiren_emri: constants OMEGA_M, H and the digitised Dalya
completeness table. Nothing under darksiren_emri/ is changed (R24, C-26).

usage: python -m analysis.rank1.budget --seed 192000 --out results/rank1_20261110
"""

import argparse
import csv
import datetime as _dt
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.interpolate import RegularGridInterpolator
from scipy.special import ndtr

from darksiren_emri.constants import OMEGA_M, H
from darksiren_emri.galaxy_catalogue.glade_completeness import (
    _COMPLETENESS_PCT,
    _DISTANCE_MPC,
)

Arr = npt.NDArray[np.float64]

C_KMS = 299792.458
SNR_THR = 20.0
TABLES = {
    "seed61000": (
        "results/campaign51_20260728/realistic_20260729/seed61000/prepared_cramer_rao_bounds.csv",
        "9a1f2a14384a9281c97ca3be312ddaab",
    ),
    "seed62000": (
        "results/campaign51_20260728/realistic_20260729/seed62000/prepared_cramer_rao_bounds.csv",
        "0dc925a59a5983afa797e1024fdec42e",
    ),
}

# Registered fiducials and priors (PREREGISTRATION_RANK1.md section 3).
PARAMS = ["h", "kappa", "s", "a", "sig_rel", "Om"]
FID = {"h": float(H), "kappa": 0.0, "s": -0.43, "a": 7.45, "sig_rel": 0.55, "Om": float(OMEGA_M)}
STEP = {"h": 0.002, "kappa": 0.02, "s": 0.02, "a": 0.02, "sig_rel": 0.02, "Om": 0.005}
PRIOR = {"kappa": 0.5, "s": 0.3, "a": 0.5, "sig_rel": 0.2, "Om": 0.02}
BETA_REL = 1.05
# Baldry et al. 2012 double Schechter (per dex, Mpc^-3).
SMF_LOGMC, SMF_PHI1, SMF_A1, SMF_PHI2, SMF_A2 = 10.66, 3.96e-3, -0.35, 0.79e-3, -1.47
LOGMS = np.linspace(7.0, 12.5, 441)
LOGM_BH = np.linspace(3.5, 7.5, 161)
Z_BINS = [0.0, 0.25, 0.5, 0.75, np.inf]
NZ_EVENT = 161
NZ_POP = 600
N_REAL = 8
SEED0 = 192000


# --------------------------------------------------------------------------- cosmology
class Cosmo:
    """Flat LCDM distances in units of c/(100 km/s/Mpc) (h-free), Mpc."""

    def __init__(self, om: float, zmax: float = 4.0, n: int = 8001) -> None:
        self.om = om
        self.zg = np.linspace(0.0, zmax, n)
        ez = np.sqrt(om * (1 + self.zg) ** 3 + 1 - om)
        dc = np.concatenate([[0.0], np.cumsum(0.5 * (1 / ez[1:] + 1 / ez[:-1]) * np.diff(self.zg))])
        self.dc = dc * C_KMS / 100.0  # Mpc * h
        self.ez = ez

    def Dc(self, z: Arr) -> Arr:
        return np.asarray(np.interp(z, self.zg, self.dc))

    def E(self, z: Arr) -> Arr:
        return np.asarray(np.sqrt(self.om * (1 + z) ** 3 + 1 - self.om))

    def dL(self, z: Arr, h: float) -> Arr:
        """Luminosity distance in Gpc."""
        return (1 + z) * self.Dc(z) / h / 1000.0

    def dVdzdOmega(self, z: Arr, h: float) -> Arr:
        """Comoving volume per unit z per steradian, Mpc^3."""
        dc = self.Dc(z) / h
        return dc**2 * (C_KMS / (100.0 * h)) / self.E(z)

    def z_of_dL(self, dl: Arr, h: float) -> Arr:
        d = (1 + self.zg) * self.dc / h / 1000.0
        return np.asarray(np.interp(dl, d, self.zg))

    def gamma(self, z: Arr) -> Arr:
        """d ln d_L / d ln(1+z)."""
        return 1 + (1 + z) * (C_KMS / 100.0) / (self.E(z) * self.Dc(z))


# --------------------------------------------------------------------------- events
@dataclass
class Events:
    name: str
    Mz: Arr
    dL: Arr
    sig_dL: Arr
    snr: Arr
    aeff: Arr  # steradian
    z: Arr = field(default_factory=lambda: np.zeros(0))

    @property
    def n(self) -> int:
        return int(len(self.dL))


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def load_events(repo: Path, name: str, cosmo: Cosmo) -> Events:
    rel, pin = TABLES[name]
    path = repo / rel
    got = md5(path)
    if got != pin:
        raise SystemExit(f"STOP: md5 mismatch for {rel}: {got} != {pin} (C-15)")
    rows = list(csv.DictReader(path.open()))

    def col(k: str) -> Arr:
        return np.array([float(r[k]) for r in rows])

    th = col("qS")
    det = col("delta_qS_delta_qS") * col("delta_phiS_delta_phiS") - col("delta_phiS_delta_qS") ** 2
    aeff = 4 * np.pi * np.abs(np.sin(th)) * np.sqrt(np.clip(det, 0, None))
    ev = Events(
        name=name,
        Mz=col("M"),
        dL=col("luminosity_distance"),
        sig_dL=np.sqrt(col("delta_luminosity_distance_delta_luminosity_distance")),
        snr=col("SNR"),
        aeff=aeff,
    )
    ev.z = cosmo.z_of_dL(ev.dL, FID["h"])
    return ev


# --------------------------------------------------------------------------- population
def phi_dex(logm: Arr, s: float) -> Arr:
    """EMRI-weighted MBH mass function per dex (unnormalised), kink at 1e5, support [1e4,1e7]."""
    out = np.where(logm >= 5.0, 10 ** (s * (logm - 5.0)), 10 ** ((s + 0.5) * (logm - 5.0)))
    return np.where((logm >= 4.0) & (logm <= 7.0), out, 0.0)


def phi_norm(s: float) -> float:
    lm = np.linspace(4.0, 7.0, 601)
    return float(np.trapezoid(phi_dex(lm, s), lm))


def smf(logms: Arr) -> Arr:
    x = 10 ** (logms - SMF_LOGMC)
    return np.asarray(
        np.log(10) * np.exp(-x) * (SMF_PHI1 * x ** (SMF_A1 + 1) + SMF_PHI2 * x ** (SMF_A2 + 1))
    )


def fcat_table(a: float, sig: float) -> RegularGridInterpolator:
    """f_cat(log M_BH, log M*_lim) = P(log M* > lim | M_BH) under SMF x relation."""
    lw = (
        np.log(smf(LOGMS))[None, :]
        - 0.5 * ((LOGM_BH[:, None] - (a + BETA_REL * (LOGMS[None, :] - 11.0))) / sig) ** 2
    )
    w = np.exp(lw - lw.max(axis=1, keepdims=True))
    cum = np.cumsum(w[:, ::-1], axis=1)[:, ::-1]
    frac = cum / cum[:, :1]
    return RegularGridInterpolator((LOGM_BH, LOGMS), frac, bounds_error=False, fill_value=None)


def fcat(tab: RegularGridInterpolator, logm: Arr, loglim: Arr) -> Arr:
    lm = np.clip(logm, LOGM_BH[0], LOGM_BH[-1])
    ll = np.clip(loglim, LOGMS[0], LOGMS[-1])
    pts = np.stack([lm.ravel(), ll.ravel()], axis=-1)
    return np.asarray(np.clip(tab(pts).reshape(lm.shape), 0.0, 1.0))


def n_cat(loglim: Arr) -> Arr:
    """Number density of catalogued galaxies above the mass limit, Mpc^-3."""
    cum = np.concatenate(
        [[0.0], np.cumsum(0.5 * (smf(LOGMS[1:]) + smf(LOGMS[:-1])) * np.diff(LOGMS))]
    )
    tot = cum[-1]
    return np.asarray(tot - np.interp(np.clip(loglim, LOGMS[0], LOGMS[-1]), LOGMS, cum))


# --------------------------------------------------------------------------- selection
@dataclass
class Selection:
    c0: float
    c1: float
    c2: float
    sig: float
    mass_dependent: bool = True

    def mu(self, lnmz: Arr) -> Arr:
        if not self.mass_dependent:
            return np.full_like(lnmz, self.c0)
        x = lnmz - np.log(5e5)
        return np.asarray(self.c0 + self.c1 * x + self.c2 * x**2)

    def pdet(self, dl: Arr, lnmz: Arr) -> Arr:
        return np.asarray(ndtr((self.mu(lnmz) - np.log(dl)) / self.sig))


def ppop(z: Arr, kappa_tot: float, cosmo: Cosmo, h: float) -> Arr:
    return cosmo.dVdzdOmega(z, h) * (1 + z) ** (kappa_tot - 1.0)


def fit_selection(ev: Events, cosmo: Cosmo, kappa0: float) -> Selection:
    lnmz = np.log(ev.Mz)
    lnhor = np.log(ev.dL * ev.snr / SNR_THR)
    zh = cosmo.z_of_dL(np.exp(lnhor), FID["h"])
    zg = np.linspace(1e-4, 4.0, 4000)
    cum = np.cumsum(ppop(zg, kappa0, cosmo, FID["h"])) * (zg[1] - zg[0])
    w = 1.0 / np.interp(zh, zg, cum)
    x = lnmz - np.log(5e5)
    A = np.stack([np.ones_like(x), x, x**2], axis=1)
    W = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(A * W, lnhor * W[:, 0], rcond=None)
    res = lnhor - A @ coef
    sig = float(np.sqrt(np.sum(w * res**2) / np.sum(w)))
    return Selection(float(coef[0]), float(coef[1]), float(coef[2]), sig)


def alpha(par: dict[str, float], sel: Selection, cosmo: Cosmo, kappa0: float, zmax: float) -> float:
    zg = np.linspace(1e-3, zmax, NZ_POP)
    lm = np.linspace(4.0, 7.0, 121)
    pz = ppop(zg, kappa0 + par["kappa"], cosmo, par["h"])
    ph = phi_dex(lm, par["s"]) / phi_norm(par["s"])
    dl = cosmo.dL(zg, par["h"])
    lnmz = np.log(10 ** lm[None, :] * (1 + zg[:, None]))
    pd = sel.pdet(dl[:, None] * np.ones_like(lnmz), lnmz)
    inner = np.trapezoid(ph[None, :] * pd, lm, axis=1)
    return float(np.trapezoid(pz * inner, zg))


def fit_kappa0(ev: Events, cosmo: Cosmo, zmax: float, n_iter: int = 2) -> tuple[float, Selection]:
    kappa0 = 0.0
    sel = fit_selection(ev, cosmo, kappa0)
    grid = np.linspace(-3.0, 5.0, 161)
    for _ in range(n_iter):
        ll = []
        base = {**FID, "kappa": 0.0}
        for k in grid:
            lp = np.log(ppop(ev.z, k, cosmo, FID["h"])) + np.log(
                sel.pdet(ev.dL, np.log(ev.Mz)) + 1e-300
            )
            ll.append(lp.sum() - ev.n * np.log(alpha(base, sel, cosmo, k, zmax)))
        kappa0 = float(grid[int(np.argmax(ll))])
        sel = fit_selection(ev, cosmo, kappa0)
    return kappa0, sel


# --------------------------------------------------------------------------- catalogue scenarios
def dalya_mass_limit(cosmo: Cosmo) -> tuple[Arr, Arr]:
    """log M*_lim(z) for GLADE+: mass-weighted completeness above the limit = Dalya f(d_L)."""
    zg = np.linspace(1e-4, 4.0, 2000)
    dl_mpc = cosmo.dL(zg, FID["h"]) * 1000.0
    f = np.interp(dl_mpc, np.array(_DISTANCE_MPC), np.array(_COMPLETENESS_PCT) / 100.0)
    mw = 10**LOGMS * smf(LOGMS)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (mw[1:] + mw[:-1]) * np.diff(LOGMS))])
    above = (cum[-1] - cum) / cum[-1]  # decreasing in LOGMS
    f = np.clip(f, 1e-9, 1.0)
    loglim = np.interp(-f, -above, LOGMS)
    loglim = np.where(f >= 0.999, 7.0, loglim)
    return zg, loglim


@dataclass
class Scenario:
    name: str
    sigz0: float  # sigma_z = sigz0 (1+z)
    outlier: float
    n_spec: int
    zg_lim: Arr
    loglim: Arr

    def lim(self, z: Arr) -> Arr:
        return np.asarray(np.interp(z, self.zg_lim, self.loglim))


def make_scenarios(cosmo: Cosmo) -> list[Scenario]:
    zg, gl = dalya_mass_limit(cosmo)
    dl05 = float(cosmo.dL(np.array([0.5]), FID["h"])[0])
    ll = np.maximum(8.5 + 2 * np.log10(cosmo.dL(zg, FID["h"]) / dl05), 7.0)
    return [
        Scenario("GLADE+ today", 0.033, 0.0, 0, zg, gl),
        Scenario("LSST+Euclid photometry", 0.04, 0.05, 0, zg, ll),
        Scenario("pointed spectroscopy (N=20)", 0.04, 0.05, 20, zg, ll),
    ]


# --------------------------------------------------------------------------- per-event grids
@dataclass
class Grids:
    zg: Arr  # (N, NZ)
    comb: dict[str, Arr]  # scenario -> C_i(z) (N, NZ)
    is_spec: dict[str, npt.NDArray[np.bool_]]
    dhat: Arr
    census: dict[str, dict[str, float]]


def event_zgrid(ev: Events, cosmo: Cosmo) -> Arr:
    gam = cosmo.gamma(ev.z)
    sz = (ev.sig_dL / ev.dL) * (1 + ev.z) / gam
    # +-10 sigma_z of the GW likelihood plus 1 % of z (covers the finite-difference h steps);
    # relative, so that the sharpest events (sigma_dL/d_L ~ 1e-3) stay resolved (dz <= sigma_z/4)
    half = 10.0 * sz + 0.01 * ev.z
    lo = np.maximum(ev.z - half, 1e-4)
    hi = ev.z + half
    t = np.linspace(0.0, 1.0, NZ_EVENT)
    return lo[:, None] + (hi - lo)[:, None] * t[None, :]


def build_comb(
    ev: Events, zg: Arr, sc: Scenario, cosmo: Cosmo, rng: np.random.Generator, zmax: float
) -> tuple[Arr, npt.NDArray[np.bool_], dict[str, float]]:
    """Catalogue contrast C_i(z) = comb / E[comb] on each event's z grid (theta-independent)."""
    tab = fcat_table(FID["a"], FID["sig_rel"])
    msrc = np.log10(ev.Mz / (1 + ev.z))
    host_cat = rng.random(ev.n) < fcat(tab, msrc, sc.lim(ev.z))
    C = np.ones_like(zg)
    counts = np.zeros(ev.n)
    # pointed spectroscopy: the n_spec events with the smallest expected candidate count
    sz_gw = (ev.sig_dL / ev.dL) * (1 + ev.z) / cosmo.gamma(ev.z)
    kexp = n_cat(sc.lim(ev.z)) * ev.aeff * cosmo.dVdzdOmega(ev.z, FID["h"]) * 4 * sz_gw
    spec = np.zeros(ev.n, dtype=bool)
    if sc.n_spec > 0:
        spec[np.argsort(kexp)[: sc.n_spec]] = True
    # outlier background: total catalogued count in the area over [0, zmax], spread uniformly
    zall = np.linspace(1e-3, zmax, 800)
    for i in range(ev.n):
        z = zg[i]
        if spec[i]:
            sp = 0.001 * (1 + ev.z[i])
            zh = ev.z[i] + sp * rng.standard_normal()
            k = np.exp(-0.5 * ((z - zh) / sp) ** 2) / (np.sqrt(2 * np.pi) * sp)
            nbar = n_cat(sc.lim(z)) * ev.aeff[i] * cosmo.dVdzdOmega(z, FID["h"])
            C[i] = k / np.maximum(nbar, 1e-300)
            counts[i] = 1
            continue
        sz = sc.sigz0 * (1 + ev.z[i])
        zlo, zhi = max(z[0] - 6 * sz, 1e-4), z[-1] + 6 * sz
        zb = np.linspace(zlo, zhi, 401)
        dz = zb[1] - zb[0]
        lam = (
            n_cat(sc.lim(zb)) * ev.aeff[i] * cosmo.dVdzdOmega(zb, FID["h"]) * dz * (1 - sc.outlier)
        )
        # photo-z scatter of a Poisson field is a Poisson field with the smoothed mean
        sb = sc.sigz0 * (1 + zb)
        Kmat = np.exp(-0.5 * ((zb[:, None] - zb[None, :]) / sb[None, :]) ** 2) / sb[None, :]
        Kmat /= Kmat.sum(axis=0, keepdims=True)
        lam_hat = Kmat @ lam
        if sc.outlier > 0:
            tot = (
                np.trapezoid(n_cat(sc.lim(zall)) * cosmo.dVdzdOmega(zall, FID["h"]), zall)
                * ev.aeff[i]
            )
            lam_hat = lam_hat + sc.outlier * tot / zmax * dz
        m = rng.poisson(lam_hat).astype(float)
        # census analogue of the production 1D ball: zhat within the z range spanned by
        # d_L +- 3 sigma over the h prior [0.60, 0.86] (C-4 grid edges)
        zb_lo = cosmo.z_of_dL(np.array([max(ev.dL[i] - 3 * ev.sig_dL[i], 1e-4)]), 0.60)[0]
        zb_hi = cosmo.z_of_dL(np.array([ev.dL[i] + 3 * ev.sig_dL[i]]), 0.86)[0]
        counts[i] = m[(zb >= zb_lo) & (zb <= zb_hi)].sum()
        zhat_host = None
        if host_cat[i]:
            if sc.outlier > 0 and rng.random() < sc.outlier:
                zhat_host = rng.uniform(0, zmax)
            else:
                zhat_host = ev.z[i] + sz * rng.standard_normal()
        # comb(z) = sum_j m_j N(z; zhat_j, sigma(zhat_j)) ; expectation with lam_hat
        sj = sc.sigz0 * (1 + zb)
        G = np.exp(-0.5 * ((z[:, None] - zb[None, :]) / sj[None, :]) ** 2) / (
            np.sqrt(2 * np.pi) * sj[None, :]
        )
        comb = G @ m
        expc = G @ lam_hat
        if zhat_host is not None:
            sh = sc.sigz0 * (1 + zhat_host)
            comb = comb + np.exp(-0.5 * ((z - zhat_host) / sh) ** 2) / (np.sqrt(2 * np.pi) * sh)
            counts[i] += float(zb_lo <= zhat_host <= zb_hi)
        ok = expc > 1e-12
        C[i] = np.where(ok, comb / np.where(ok, expc, 1.0), 1.0)
    nonempty = counts[counts > 0]
    census = {
        "frac_empty": float(np.mean(counts == 0)),
        "n_empty": float(np.sum(counts == 0)),
        "median_K_nonempty": float(np.median(nonempty)) if len(nonempty) else 0.0,
        "mean_K_nonempty": float(np.mean(nonempty)) if len(nonempty) else 0.0,
        "host_catalogued_share": float(np.mean(host_cat | spec)),
        "host_catalogued_n": float(np.sum(host_cat | spec)),
    }
    return C, spec, census


# --------------------------------------------------------------------------- likelihood
@dataclass
class Model:
    ev: Events
    sel: Selection
    kappa0: float
    zmax: float
    scen: Scenario | None  # None => population only (C = 1)
    mass_kink: bool = True


def per_event_lnnum(
    par: dict[str, float],
    m: Model,
    zg: Arr,
    C: Arr | None,
    spec: npt.NDArray[np.bool_] | None,
    dhat: Arr,
) -> Arr:
    cosmo = Cosmo(par["Om"]) if abs(par["Om"] - FID["Om"]) > 1e-12 else _COSMO_FID
    ev = m.ev
    pz = ppop(zg, m.kappa0 + par["kappa"], cosmo, par["h"])
    logm = np.log10(ev.Mz)[:, None] - np.log10(1 + zg)
    s = par["s"]
    if m.mass_kink:
        ph = phi_dex(logm, s) / phi_norm(s)
    else:
        ph = 10 ** (s * (logm - 5.0)) / _pl_norm(s)
    rho = pz * ph
    if C is not None and m.scen is not None:
        tab = fcat_table(par["a"], par["sig_rel"])
        f = fcat(tab, logm, m.scen.lim(zg))
        if spec is not None:
            f = np.where(spec[:, None], 1.0, f)
        rho = rho * (f * C + 1.0 - f)
    dl = cosmo.dL(zg, par["h"])
    gw = np.exp(-0.5 * ((dhat[:, None] - dl) / ev.sig_dL[:, None]) ** 2)
    pd = m.sel.pdet(dl, np.log(ev.Mz)[:, None] * np.ones_like(dl))
    num = np.trapezoid(rho * gw * pd, zg, axis=1)
    return np.log(np.maximum(num, 1e-300))


def _pl_norm(s: float) -> float:
    lm = np.linspace(4.0, 7.0, 601)
    return float(np.trapezoid(10 ** (s * (lm - 5.0)), lm))


def lnalpha(par: dict[str, float], m: Model) -> float:
    cosmo = Cosmo(par["Om"]) if abs(par["Om"] - FID["Om"]) > 1e-12 else _COSMO_FID
    if m.mass_kink:
        return float(np.log(alpha(par, m.sel, cosmo, m.kappa0, m.zmax)))
    zg = np.linspace(1e-3, m.zmax, NZ_POP)
    lm = np.linspace(4.0, 7.0, 121)
    pz = ppop(zg, m.kappa0 + par["kappa"], cosmo, par["h"])
    ph = 10 ** (par["s"] * (lm - 5.0)) / _pl_norm(par["s"])
    dl = cosmo.dL(zg, par["h"])
    lnmz = np.log(10 ** lm[None, :] * (1 + zg[:, None]))
    pd = m.sel.pdet(dl[:, None] * np.ones_like(lnmz), lnmz)
    return float(np.log(np.trapezoid(pz * np.trapezoid(ph[None, :] * pd, lm, axis=1), zg)))


_COSMO_FID = Cosmo(FID["Om"])


def hessian_terms(
    m: Model,
    params: list[str],
    zg: Arr,
    C: Arr | None,
    spec: npt.NDArray[np.bool_] | None,
    dhat: Arr,
) -> tuple[Arr, Arr, Arr, Arr]:
    """Per-event negative Hessian of ln num_i and the negative Hessian of ln alpha.

    Returns (H_ev (N, P, P), H_alpha (P, P), G_ev (N, P), G_alpha (P)); the Fisher of any subset
    S is sum_{i in S} H_ev[i] - |S| * H_alpha. G are first derivatives (central differences), used
    by the post-hoc outer-product-of-scores diagnostic.
    """
    P = len(params)
    cache: dict[tuple[int, ...], tuple[Arr, float]] = {}

    def f(offs: tuple[int, ...]) -> tuple[Arr, float]:
        if offs not in cache:
            par = dict(FID)
            for k, o in zip(params, offs):
                par[k] = FID[k] + o * STEP[k]
            cache[offs] = (per_event_lnnum(par, m, zg, C, spec, dhat), lnalpha(par, m))
        return cache[offs]

    zero = tuple([0] * P)
    Hev = np.zeros((m.ev.n, P, P))
    Ha = np.zeros((P, P))
    Gev = np.zeros((m.ev.n, P))
    Ga = np.zeros(P)
    f0, a0 = f(zero)
    for j in range(P):
        for k in range(j, P):
            if j == k:
                up = list(zero)
                up[j] = 1
                dn = list(zero)
                dn[j] = -1
                fu, au = f(tuple(up))
                fd, ad = f(tuple(dn))
                d2 = (fu - 2 * f0 + fd) / STEP[params[j]] ** 2
                d2a = (au - 2 * a0 + ad) / STEP[params[j]] ** 2
                Gev[:, j] = (fu - fd) / (2 * STEP[params[j]])
                Ga[j] = (au - ad) / (2 * STEP[params[j]])
            else:
                vals = []
                for sj, sk in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    o = list(zero)
                    o[j] = sj
                    o[k] = sk
                    vals.append(f(tuple(o)))
                d2 = (vals[0][0] - vals[1][0] - vals[2][0] + vals[3][0]) / (
                    4 * STEP[params[j]] * STEP[params[k]]
                )
                d2a = (vals[0][1] - vals[1][1] - vals[2][1] + vals[3][1]) / (
                    4 * STEP[params[j]] * STEP[params[k]]
                )
            Hev[:, j, k] = Hev[:, k, j] = -d2
            Ha[j, k] = Ha[k, j] = -d2a
    return Hev, Ha, Gev, Ga


def fisher_opg(Gev: Arr, Ga: Arr, mask: npt.NDArray[np.bool_] | None = None) -> Arr:
    """Outer product of per-event scores (post-hoc diagnostic; PSD by construction)."""
    if mask is None:
        mask = np.ones(len(Gev), dtype=bool)
    s = Gev[mask] - Ga[None, :]
    return np.asarray(s.T @ s)


def fisher(Hev: Arr, Ha: Arr, mask: npt.NDArray[np.bool_] | None = None) -> Arr:
    if mask is None:
        mask = np.ones(len(Hev), dtype=bool)
    return np.asarray(Hev[mask].sum(axis=0) - mask.sum() * Ha)


def sigmas(
    F: Arr, params: list[str], prior_scale: float = 1.0, prior_over: dict[str, float] | None = None
) -> dict[str, float]:
    """sigma_h conditional, marginalised, ratio."""
    ih = params.index("h")
    cond = 1.0 / np.sqrt(F[ih, ih]) if F[ih, ih] > 0 else float("nan")
    Pm = np.zeros_like(F)
    for j, k in enumerate(params):
        if k == "h":
            continue
        w = (prior_over or {}).get(k, PRIOR[k]) * prior_scale
        Pm[j, j] = 1.0 / w**2
    M = F + Pm
    try:
        ev = np.linalg.eigvalsh(M)
        if ev.min() <= 0:
            return {
                "cond": float(cond),
                "marg": float("nan"),
                "R": float("nan"),
                "posdef": 0.0,
                "min_eig": float(ev.min()),
                "max_eig": float(ev.max()),
            }
        marg = float(np.sqrt(np.linalg.inv(M)[ih, ih]))
    except np.linalg.LinAlgError:
        return {"cond": float(cond), "marg": float("nan"), "R": float("nan"), "posdef": 0.0}
    return {"cond": float(cond), "marg": marg, "R": marg / cond, "posdef": 1.0}


# --------------------------------------------------------------------------- driver
def realisation(
    ev: Events,
    cosmo: Cosmo,
    scens: list[Scenario],
    sel: Selection,
    kappa0: float,
    zmax: float,
    seed: int,
    params: list[str],
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    zg = event_zgrid(ev, cosmo)
    dhat = ev.dL + ev.sig_dL * rng.standard_normal(ev.n)
    out: dict[str, Any] = {"census": {}}
    pop = Model(ev, sel, kappa0, zmax, None)
    out["POP"] = hessian_terms(pop, params, zg, None, None, dhat)
    for sc in scens:
        C, spec, census = build_comb(ev, zg, sc, cosmo, rng, zmax)
        out["census"][sc.name] = census
        mm = Model(ev, sel, kappa0, zmax, sc)
        out[sc.name] = hessian_terms(mm, params, zg, C, spec, dhat)
    return out


def summarise(
    reals: list[dict[str, Any]],
    keys: list[str],
    params: list[str],
    masks: dict[str, Any],
    prior_scale: float = 1.0,
    prior_over: dict[str, float] | None = None,
    n_boot: int = 200,
) -> dict[str, Any]:
    res: dict[str, Any] = {}
    rngb = np.random.default_rng(SEED0 + 999)
    for key in keys:
        res[key] = {}
        for mname, mask in masks.items():
            Fs = [fisher(r[key][0], r[key][1], mask) for r in reals]
            Fp = [fisher(r["POP"][0], r["POP"][1], mask) for r in reals]
            F = np.mean(Fs, axis=0)
            FP = np.mean(Fp, axis=0)
            full = sigmas(F, params, prior_scale, prior_over)
            popu = sigmas(FP, params, prior_scale, prior_over)
            boots = []
            for _ in range(n_boot):
                idx = rngb.integers(0, len(reals), len(reals))
                boots.append(
                    sigmas(np.mean([Fs[i] for i in idx], axis=0), params, prior_scale, prior_over)[
                        "R"
                    ]
                )
            se = float(np.nanstd(boots))

            def assoc(a: float, b: float) -> float:
                d = 1 / a**2 - 1 / b**2
                return float(1 / np.sqrt(d)) if d > 0 else float("nan")

            Fo = np.mean([fisher_opg(r[key][2], r[key][3], mask) for r in reals], axis=0)
            Fpo = np.mean([fisher_opg(r["POP"][2], r["POP"][3], mask) for r in reals], axis=0)
            Fn = F + np.diag(
                [
                    0.0
                    if k == "h"
                    else 1.0 / ((prior_over or {}).get(k, PRIOR[k]) * prior_scale) ** 2
                    for k in params
                ]
            )
            ew, evec = np.linalg.eigh(Fn)
            res[key][mname + "__diag"] = {
                "FULL_OPG": sigmas(Fo, params, prior_scale, prior_over),
                "POPULATION_OPG": sigmas(Fpo, params, prior_scale, prior_over),
                "eig_F_plus_P": [float(x) for x in ew],
                "eigvec_min": {k: float(v) for k, v in zip(params, evec[:, 0])},
                "F_diag": {k: float(F[i, i]) for i, k in enumerate(params)},
            }
            res[key][mname] = {
                "FULL": {**full, "R_se": se},
                "POPULATION": popu,
                "ASSOCIATION": {
                    "cond": assoc(full["cond"], popu["cond"]),
                    "marg": assoc(full["marg"], popu["marg"]),
                },
                "n_events": int(np.sum(mask)) if mask is not None else int(len(reals[0]["POP"][0])),
                "pop_share_cond_info": float((1 / popu["cond"] ** 2) / (1 / full["cond"] ** 2)),
            }
    return res


def git_commit(repo: Path) -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True
        ).stdout.strip()
    except OSError:
        return "unknown"


def run_table(repo: Path, name: str, params: list[str], n_real: int) -> dict[str, Any]:
    cosmo = _COSMO_FID
    ev = load_events(repo, name, cosmo)
    zmax = float(1.3 * ev.z.max())
    kappa0, sel = fit_kappa0(ev, cosmo, zmax)
    scens = make_scenarios(cosmo)
    reals = [
        realisation(ev, cosmo, scens, sel, kappa0, zmax, SEED0 + r, params) for r in range(n_real)
    ]
    return {"ev": ev, "reals": reals, "kappa0": kappa0, "sel": sel, "scens": scens, "zmax": zmax}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=SEED0)
    ap.add_argument("--out", type=Path, default=Path("results/rank1_20261110"))
    ap.add_argument("--n_real", type=int, default=N_REAL)
    ap.add_argument("--repo", type=Path, default=Path("."))
    args = ap.parse_args(argv)
    if args.seed != SEED0:
        raise SystemExit("the registered seed base is 192000 (PREREGISTRATION_RANK1.md section 1)")
    out: Path = args.out
    out.mkdir(parents=True, exist_ok=True)
    meta = {
        "git_commit": git_commit(args.repo),
        "timestamp": _dt.datetime.now(_dt.UTC).isoformat(),
        "argv": sys.argv,
        "seed_base": SEED0,
        "n_real": args.n_real,
        "inputs": {k: v[1] for k, v in TABLES.items()},
        "prereg": ".planning/missions/M192/PREREGISTRATION_RANK1.md",
    }
    (out / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    P5 = PARAMS[:5]
    result: dict[str, Any] = {"meta": meta}
    keys = ["GLADE+ today", "LSST+Euclid photometry", "pointed spectroscopy (N=20)"]
    for tname in ("seed61000", "seed62000"):
        print(f"[{tname}] building realisations", flush=True)
        T = run_table(args.repo, tname, P5, args.n_real)
        ev = T["ev"]
        masks: dict[str, Any] = {"all": np.ones(ev.n, dtype=bool)}
        for lo, hi in zip(Z_BINS[:-1], Z_BINS[1:]):
            masks[f"z[{lo},{hi})"] = (ev.z >= lo) & (ev.z < hi)
        block: dict[str, Any] = {
            "kappa0": T["kappa0"],
            "selection": T["sel"].__dict__,
            "zmax": T["zmax"],
            "census": T["reals"][0]["census"],
            "census_mean": {
                k: {
                    kk: float(np.mean([r["census"][k][kk] for r in T["reals"]]))
                    for kk in T["reals"][0]["census"][k]
                }
                for k in keys
            },
            "registered": summarise(T["reals"], keys, P5, masks),
        }
        if tname == "seed61000":
            sens: dict[str, Any] = {}
            sens["priors_x0.5"] = summarise(
                T["reals"], keys, P5, {"all": masks["all"]}, prior_scale=0.5
            )
            sens["priors_x2"] = summarise(
                T["reals"], keys, P5, {"all": masks["all"]}, prior_scale=2.0
            )
            sens["kappa_prior_0.25"] = summarise(
                T["reals"], keys, P5, {"all": masks["all"]}, prior_over={"kappa": 0.25}
            )
            sens["kappa_prior_1.0"] = summarise(
                T["reals"], keys, P5, {"all": masks["all"]}, prior_over={"kappa": 1.0}
            )
            # Pozzoli-like count: N = 139, median over 10 random subsets
            subs = []
            for j in range(10):
                rj = np.random.default_rng(SEED0 + 100 + j)
                msk = np.zeros(ev.n, dtype=bool)
                msk[rj.choice(ev.n, 139, replace=False)] = True
                subs.append(summarise(T["reals"], keys, P5, {"N139": msk}, n_boot=20))
            poz: dict[str, Any] = {}
            for k in keys:
                poz[k] = {
                    ch: {
                        q: float(np.nanmedian([s[k]["N139"][ch][q] for s in subs]))
                        for q in ("cond", "marg") + (("R",) if ch != "ASSOCIATION" else ())
                    }
                    for ch in ("FULL", "POPULATION", "ASSOCIATION")
                }
            sens["pozzoli_N139"] = poz
            block["sensitivity"] = sens
            # Omega_m free (sensitivity row): separate 6-parameter Hessians, 2 realisations
            print("[seed61000] Omega_m-free row", flush=True)
            reals6 = [
                realisation(
                    ev, _COSMO_FID, T["scens"], T["sel"], T["kappa0"], T["zmax"], SEED0 + r, PARAMS
                )
                for r in range(2)
            ]
            block["sensitivity"]["Om_free"] = summarise(
                reals6, keys, PARAMS, {"all": masks["all"]}, n_boot=20
            )
        result[tname] = block
        (out / "budget.json").write_text(json.dumps(result, indent=2, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
