# PROVENANCE (mission M192, 2026-11-10): verbatim copy of the outside advisor's toy
# source: garden wiki/assets/darksiren-advisor-2026-11-10/toy/q4_line_toy.py
# garden commit 55a3be19dd4cbcfe0ca4a48b92b0997d197de72a, source md5 d3e8c73c0a3d09f502e27471137189c4
# (lines below this header are byte-identical to the source; the mass-dimension extension is
# analysis/rank1/toy/q4_mass_toy.py). Not run by M192; kept as the reference implementation.
# ruff: noqa
"""Q4 toy: is a catalogue + completion dark-siren likelihood calibrated at ~5 % in-catalogue share?

Sky = P independent patches (each = one LISA error box, a 1-D line in z). Galaxies Poisson in each
patch with intensity nb(z) = NBOX dVc/dz / int (flat LCDM, Om=0.2726, z<0.6); Schechter
luminosities (alpha=-1.07, L>0.01L*). Catalogue = flux-limited at h_true: L/dL^2 > F.
Hosts: number-weighted over all galaxies of the sky (every galaxy equally likely), kept if detected:
deterministic true dL < DH (= dL(z=0.4)). Data dhat = dL_true(1+0.04 eps). Exact (spec) catalogue z.
Gray-type likelihood, number-weighted, h-dependent completeness, theta-deterministic selection:
  num_i(h) = sum_{g in cat, patch i} N(dhat; dL(z_g,h)) 1[dL<DH] + int dz nb Punc(z,h) N(dhat;dL(z,h)) 1[dL<DH]
  alpha(h) = (1/P) #{g in sky catalogue: dL(z_g,h)<DH} + int dz nb Punc(z,h) 1[dL(z,h)<DH]
  Punc(z,h) = P_number(L < F dL(z,h)^2)
Variants: COMPLETE (every galaxy catalogued; same events), LUMW (estimator uses luminosity-weighted
completeness for the number-weighted hosts), FROZEN (Punc evaluated at h_fid = 0.70).
usage: python3 q4_line_toy.py SEED FLIM NUNIVERSES
"""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
FLIM = float(sys.argv[2]) if len(sys.argv) > 2 else 2e-6
NU = int(sys.argv[3]) if len(sys.argv) > 3 else 40
OM, C, HT, SIG, NEV, P, NBOX = 0.2726, 299792.458, 0.73, 0.04, 300, 600, 1000.0
zg = np.linspace(0, 0.6, 3001); Ez = np.sqrt(OM * (1 + zg) ** 3 + 1 - OM)
Dc = np.concatenate([[0], np.cumsum(0.5 * (1 / Ez[1:] + 1 / Ez[:-1]) * np.diff(zg))]) * C / 100
zg, Ez, Dc = zg[1:], Ez[1:], Dc[1:]
dvdz = Dc ** 2 / Ez
nb = NBOX * dvdz / np.trapezoid(dvdz, zg)
def dLz(z, h): return (1 + z) * np.interp(z, zg, Dc) / h
DH = dLz(0.40, HT)
lx = np.linspace(-2, 1.5, 3501); xx = 10 ** lx
dn = xx ** (-1.07 + 1) * np.exp(-xx)
cumN = np.cumsum(dn[::-1])[::-1]; cumN /= cumN[0]
cumL = np.cumsum((xx * dn)[::-1])[::-1]; cumL /= cumL[0]
def p_unc(xlim, lum=False):
    return 1 - np.interp(np.log10(np.clip(xlim, 1e-2, 10 ** 1.5)), lx, cumL if lum else cumN)
hgrid = np.linspace(0.43, 1.03, 241)  # symmetric about h_true: avoids prior-truncation pull

def sky():
    n = rng.poisson(NBOX, P); patch = np.repeat(np.arange(P), n)
    z = np.abs(rng.choice(zg, n.sum(), p=nb / nb.sum()) + rng.uniform(-1e-4, 1e-4, n.sum()))
    x = 10 ** rng.choice(lx, n.sum(), p=dn / dn.sum())
    incat = x > FLIM * dLz(z, HT) ** 2
    det = np.where(dLz(z, HT) < DH)[0]
    hosts = rng.choice(det, NEV, replace=False)
    dhat = dLz(z[hosts], HT) * (1 + SIG * rng.standard_normal(NEV))
    return z, patch, incat, hosts, dhat

def loglike(z, patch, cat, hosts, dhat, mode):
    LL = np.zeros(len(hgrid)); hp = patch[hosts]
    zc, pc = z[cat], patch[cat]
    order = np.argsort(pc); zc, pc = zc[order], pc[order]
    starts = np.searchsorted(pc, hp); ends = np.searchsorted(pc, hp, side="right")
    for k, h in enumerate(hgrid):
        dl = dLz(zg, h); sel = dl < DH
        pu = np.zeros_like(zg) if mode == "COMPLETE" else p_unc(FLIM * dLz(zg, 0.70 if mode == "FROZEN" else h) ** 2, lum=(mode == "LUMW"))
        dens = nb * pu * sel
        dlc = dLz(zc, h)
        wc = 2.3 if mode == "MIX2.3" else 1.0   # catalogue galaxies over-weighted x2.3 vs dark (cf. repo w_G 0.1215 vs 0.0523)
        alpha = wc * (dlc < DH).sum() / P + np.trapezoid(dens, zg)
        drk = np.trapezoid(dens[None, :] * np.exp(-0.5 * ((dhat[:, None] - dl[None, :]) / (SIG * dl[None, :])) ** 2) / dl[None, :], zg, axis=1)
        cat_t = np.zeros(NEV)
        for i in range(NEV):
            g = dlc[starts[i]:ends[i]]; g = g[g < DH]
            cat_t[i] = (np.exp(-0.5 * ((dhat[i] - g) / (SIG * g)) ** 2) / g).sum()
        LL[k] = np.log(wc * cat_t + drk + 1e-300).sum() - NEV * np.log(alpha)
    return LL

def summ(LL):
    p = np.exp(LL - LL.max()); p /= p.sum(); m = (p * hgrid).sum(); cdf = np.cumsum(p)
    return m, np.sqrt((p * (hgrid - m) ** 2).sum()), np.interp(HT, hgrid, cdf), p[0] + p[-1], np.interp(0.5, cdf, hgrid)

if __name__ == "__main__":
    res = {m: [] for m in ("CORRECT", "COMPLETE", "LUMW", "FROZEN", "MIX2.3")}; shares = []
    for u in range(NU):
        z, patch, incat, hosts, dhat = sky(); shares.append(incat[hosts].mean())
        for m in res:
            cat = np.ones_like(incat) if m == "COMPLETE" else incat
            res[m].append(summ(loglike(z, patch, cat, hosts, dhat, m)))
    print(f"FLIM={FLIM:g}: in-catalogue share of detected hosts {np.mean(shares):.3f}; catalogued fraction of all galaxies {incat.mean():.4f}; "
          f"events/universe {NEV}; universes {NU}; galaxies per box ~{NBOX:.0f}")
    for m, r in res.items():
        r = np.array(r); b = (r[:, 0] - HT) / r[:, 1]
        print(f"  {m:9s} <mean>-truth {np.mean(r[:,0])-HT:+.4f}  <sd> {np.mean(r[:,1]):.4f}  std(means) {np.std(r[:,0]):.4f}  "
              f"<pull> {np.mean(b):+.2f}+-{np.std(b)/np.sqrt(len(b)):.2f}  cov68 {np.mean((r[:,2]>.16)&(r[:,2]<.84)):.2f}  <PIT> {np.mean(r[:,2]):.2f}  <median-truth>/sd {np.mean((r[:,4]-HT)/r[:,1]):+.2f}")
    np.save(f"q4_res_{sys.argv[1]}.npy", np.array([res[m] for m in res]))
