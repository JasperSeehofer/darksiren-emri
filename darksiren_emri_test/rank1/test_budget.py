"""Limiting-case tests for the Rank-1 information budget (mission M192, item 2).

Each test is one limiting case of docs/derivations/rank1_information_budget.md section 4.
"""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import numpy.typing as npt
import pytest

_REPO = Path(__file__).resolve().parents[2]


def _load() -> ModuleType:
    path = _REPO / "analysis" / "rank1" / "budget.py"
    spec = importlib.util.spec_from_file_location("rank1_budget", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["rank1_budget"] = mod
    spec.loader.exec_module(mod)
    return mod


B: Any = _load()


def _events(n: int = 6, z0: float = 0.3, rel: float = 0.03) -> Any:
    cosmo = B._COSMO_FID
    z = np.linspace(z0 * 0.7, z0 * 1.3, n)
    dl = cosmo.dL(z, B.FID["h"])
    return B.Events(
        name="synthetic",
        Mz=np.full(n, 5e5) * (1 + z),
        dL=dl,
        sig_dL=rel * dl,
        snr=np.full(n, 40.0),
        aeff=np.full(n, 3e-4),
        z=z,
    )


def _no_selection() -> Any:
    return B.Selection(c0=50.0, c1=0.0, c2=0.0, sig=1.0, mass_dependent=False)


def _scenario() -> Any:
    zg = np.linspace(0, 4, 50)
    return B.Scenario("test", 0.04, 0.0, 0, zg, np.full_like(zg, 9.0))


def _host_comb(ev: Any, zg: npt.NDArray[np.float64], sz: float) -> npt.NDArray[np.float64]:
    """C_i(z) for a unique host at the true z with a Gaussian redshift kernel of width sz."""
    k = np.exp(-0.5 * ((zg - ev.z[:, None]) / sz) ** 2) / (np.sqrt(2 * np.pi) * sz)
    return np.asarray(k)


def test_unique_host_complete_catalogue_limit() -> None:
    """Limit 1: sigma_h/h = [sum (d_L/sigma_dL)^2]^(-1/2) with unique hosts and no selection."""
    ev = _events()
    zg = B.event_zgrid(ev, B._COSMO_FID)
    C = _host_comb(ev, zg, 1e-5)
    spec = np.ones(ev.n, dtype=bool)
    m = B.Model(ev, _no_selection(), 0.0, 2.0, _scenario())
    Hev, Ha, *_ = B.hessian_terms(m, ["h"], zg, C, spec, ev.dL.copy())
    F = B.fisher(Hev, Ha)
    expect = np.sum((ev.dL / B.FID["h"]) ** 2 / ev.sig_dL**2)
    assert F[0, 0] == pytest.approx(expect, rel=0.02)


def test_unique_host_photoz_limit() -> None:
    """Limit 2: sigma_lnh^2 = (sigma_dL/d_L)^2 + (gamma sigma_z/(1+z))^2."""
    ev = _events(n=1, z0=0.3, rel=0.03)
    sz = 0.005
    zg = B.event_zgrid(ev, B._COSMO_FID)
    C = _host_comb(ev, zg, sz)
    m = B.Model(ev, _no_selection(), 0.0, 2.0, _scenario())
    Hev, Ha, *_ = B.hessian_terms(m, ["h"], zg, C, np.ones(1, dtype=bool), ev.dL.copy())
    F = B.fisher(Hev, Ha)
    gam = float(B._COSMO_FID.gamma(ev.z)[0])
    var_lnh = 0.03**2 + (gam * sz / (1 + ev.z[0])) ** 2
    assert 1.0 / np.sqrt(F[0, 0]) / B.FID["h"] == pytest.approx(np.sqrt(var_lnh), rel=0.05)


def test_gamma_value() -> None:
    """gamma(0.47) = 3.76 at Omega_m = 0.2726 (advisor correction (b))."""
    assert float(B._COSMO_FID.gamma(np.array([0.47]))[0]) == pytest.approx(3.76, abs=0.02)


def test_empty_catalogue_equals_population() -> None:
    """Limit 3: C = 1 gives the population-only Fisher and no information on a, sigma_rel."""
    ev = _events()
    zg = B.event_zgrid(ev, B._COSMO_FID)
    sel = B.Selection(c0=np.log(3.0), c1=0.3, c2=-0.05, sig=0.4)
    params = ["h", "kappa", "a", "sig_rel"]
    full = B.Model(ev, sel, 0.0, 2.0, _scenario())
    pop = B.Model(ev, sel, 0.0, 2.0, None)
    H1, A1, *_ = B.hessian_terms(full, params, zg, np.ones_like(zg), None, ev.dL.copy())
    H0, A0, *_ = B.hessian_terms(pop, params, zg, None, None, ev.dL.copy())
    np.testing.assert_allclose(B.fisher(H1, A1), B.fisher(H0, A0), rtol=1e-9, atol=1e-6)
    F = B.fisher(H1, A1)
    assert np.all(np.abs(F[2:, :]) < 1e-6)


def test_zero_prior_width_gives_conditional() -> None:
    """Limit 4: sigma_marg -> sigma_cond as every prior width -> 0."""
    rng = np.random.default_rng(0)
    A = rng.standard_normal((5, 5))
    F = A @ A.T + 5 * np.eye(5)
    s = B.sigmas(F, B.PARAMS[:5], prior_scale=1e-6)
    assert s["marg"] == pytest.approx(s["cond"], rel=1e-6)
    s1 = B.sigmas(F, B.PARAMS[:5])
    assert s1["R"] >= 1.0


def test_slope_evolution_identity() -> None:
    """Limit 5: pure power law, mass-independent selection: F_hs = -F_hk + N d_h d_k ln alpha."""
    ev = _events()
    zg = B.event_zgrid(ev, B._COSMO_FID)
    sel = B.Selection(c0=np.log(1.5), c1=0.0, c2=0.0, sig=0.4, mass_dependent=False)
    m = B.Model(ev, sel, 0.0, 2.0, None, mass_kink=False)
    rng = np.random.default_rng(1)
    dhat = ev.dL + ev.sig_dL * rng.standard_normal(ev.n)
    Hev, Ha, *_ = B.hessian_terms(m, ["h", "kappa", "s"], zg, None, None, dhat)
    F = B.fisher(Hev, Ha)
    lhs = F[0, 2]
    rhs = -F[0, 1] + ev.n * (-Ha[0, 1])
    assert lhs == pytest.approx(rhs, rel=1e-3, abs=1e-3 * abs(F[0, 1]))


def test_input_md5_pin_stops(monkeypatch: pytest.MonkeyPatch) -> None:
    """C-15: a changed input table stops the run."""
    rel, _ = B.TABLES["seed61000"]
    monkeypatch.setitem(B.TABLES, "seed61000", (rel, "0" * 32))
    with pytest.raises(SystemExit):
        B.load_events(_REPO, "seed61000", B._COSMO_FID)


def test_comb_is_deterministic_and_unbiased() -> None:
    """Same seed -> same catalogue contrast; mean contrast over impostor fields is ~1."""
    ev = _events(n=4, z0=0.3)
    zg = B.event_zgrid(ev, B._COSMO_FID)
    sc = B.Scenario("t", 0.04, 0.0, 0, np.linspace(0, 4, 50), np.full(50, 12.4))
    c1, _, _ = B.build_comb(ev, zg, sc, B._COSMO_FID, np.random.default_rng(3), 2.0)
    c2, _, _ = B.build_comb(ev, zg, sc, B._COSMO_FID, np.random.default_rng(3), 2.0)
    np.testing.assert_array_equal(c1, c2)
    sc2 = B.Scenario("t2", 0.04, 0.0, 0, np.linspace(0, 4, 50), np.full(50, 9.0))
    ev2 = B.Events("x", ev.Mz, ev.dL, ev.sig_dL, ev.snr, np.full(ev.n, 3e-6), ev.z)
    cs = [
        B.build_comb(ev2, zg, sc2, B._COSMO_FID, np.random.default_rng(10 + k), 2.0)[0]
        for k in range(40)
    ]
    assert float(np.mean(cs)) == pytest.approx(1.0, abs=0.15)
