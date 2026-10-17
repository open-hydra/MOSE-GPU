#!/usr/bin/env python3
"""gpu-port-mcdiff unit test, independent check: mixture-averaged D_k (Curtiss-Hirschfelder, as co_DS_expr) in numpy
from diffusion.dat + phase.txt, against the host (Dh) and device (Dd) results of test_co_DS_dev.
usage: check_py.py <input dir> states.bin out.bin"""
import re, sys, numpy as np
d, fs, fo = sys.argv[1:4]
W = np.array([float(l.split()[1]) for l in open(f"{d}/phase.txt").read().splitlines()[1:] if l.strip()])
txt = open(f"{d}/diffusion.dat").read(); pref = float(re.search(r"Pref=([0-9.Ee+-]+)", txt).group(1))
Z = [np.array(z.split("\n", 2)[2].split(), float).reshape(-1, 2) for z in re.split(r'ZONE\s+T="', txt)[1:]]
assert len(Z) == 21 and all(z[0, 0] == 1.0 for z in Z)
tmax = int(Z[0][-1, 0]); tab = np.array([np.concatenate([z[:1, 1], z[:, 1]]) for z in Z])   # index 0 = copy of T = 1
with open(fs, "rb") as f:
    n = int(np.fromfile(f, "<i4", 1)[0]); R = np.fromfile(f, "<f8", 7 * n).reshape(n, 7).T
    T = np.fromfile(f, "<f8", n); P = np.fromfile(f, "<f8", n)
with open(fo, "rb") as f:
    assert int(np.fromfile(f, "<i4", 1)[0]) == n
    Dh = np.fromfile(f, "<f8", 7 * n).reshape(n, 7).T; Dd = np.fromfile(f, "<f8", 7 * n).reshape(n, 7).T
rho = R.sum(0); Wm = 1.0 / (R / rho / W[:, None]).sum(0); X = R * Wm / (rho * W[:, None])
Tl = np.floor(T).astype(int); Td = T - Tl
i1, i2 = np.clip(Tl, 0, tmax), np.clip(Tl + 1, 0, tmax)
den = np.zeros((7, n)); dsum = np.zeros((7, n)); q = 0
for i in range(7):
    for j in range(i + 1, 7):
        d0, d1 = tab[q][i1], tab[q][i2]; Dv = d0 + (d1 - d0) * Td
        den[i] += X[j] / Dv; den[j] += X[i] / Dv; dsum[i] += Dv; dsum[j] += Dv; q += 1
with np.errstate(divide="ignore", invalid="ignore"):
    Dp = np.where(1.0 - X < 1e-10, dsum / 6.0, (1.0 - X) / den) * (pref / P)
def rel(a, b): return np.abs(a - b) / np.maximum(np.abs(b), np.finfo(float).tiny)
sp = ["O2", "H2O", "H2", "H", "O", "OH", "N2"]
for lab, D in (("host", Dh), ("device", Dd)):
    r = rel(D, Dp); print(f"python vs {lab}: max rel {np.nanmax(r):.3e}  (per species " + ", ".join(f"{sp[k]} {np.nanmax(r[k]):.1e}" for k in range(7)) + f"); NaN python {np.isnan(Dp).sum()}")
print(f"host vs device: identical {np.sum(Dh == Dd)} / {Dh.size}, max rel {np.max(rel(Dd, Dh)):.3e}")
nr = n - 10016
print(f"range of D_k on the real states [m2/s]: " + ", ".join(f"{sp[k]} {Dh[k, :nr].min():.2e}-{Dh[k, :nr].max():.2e}" for k in range(7)))
print("CHECK_PY_DONE")
