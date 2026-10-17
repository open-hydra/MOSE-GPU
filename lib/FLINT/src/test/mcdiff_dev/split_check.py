#!/usr/bin/env python3
"""gpu-port-mcdiff unit test: where the host/device differences sit. Splits real (pocket) states from edge/random ones and normalises
the relative difference by the conditioning of the numerator (1 - X_k): an error of one ulp in X_k becomes eps/(1 - X_k) in D_k.
usage: split_check.py <input dir> states.bin out.bin [nreal]"""
import sys, numpy as np
d, fs, fo = sys.argv[1:4]; nreal = int(sys.argv[4]) if len(sys.argv) > 4 else 1224000
W = np.array([float(l.split()[1]) for l in open(f"{d}/phase.txt").read().splitlines()[1:] if l.strip()])
with open(fs, "rb") as f:
    n = int(np.fromfile(f, "<i4", 1)[0]); R = np.fromfile(f, "<f8", 7 * n).reshape(n, 7).T
with open(fo, "rb") as f:
    np.fromfile(f, "<i4", 1); Dh = np.fromfile(f, "<f8", 7 * n).reshape(n, 7).T; Dd = np.fromfile(f, "<f8", 7 * n).reshape(n, 7).T
rho = R.sum(0); X = (R / W[:, None]) / (R / W[:, None]).sum(0)
rel = np.abs(Dd - Dh) / np.abs(Dh); eps = np.finfo(float).eps
cond = rel * (1 - X) / eps          # difference in units of eps/(1 - X_k)
for lab, sl in (("real", slice(0, nreal)), ("edge+random", slice(nreal, n))):
    r, c, x = rel[:, sl], cond[:, sl], X[:, sl]
    print(f"{lab:12s}: non-identical {np.sum(r > 0)} / {r.size}, max rel {r.max():.2e}, max rel*(1-X)/eps {c.max():.1f}, "
          f"X_k at the max rel {x.ravel()[r.argmax()]:.6f}")
    for s, sp in enumerate(["O2", "H2O", "H2", "H", "O", "OH", "N2"]):
        k = r[s].argmax(); print(f"   {sp:4s} max rel {r[s, k]:.2e} at X {x[s, k]:.6f}; max rel*(1-X)/eps {c[s].max():.1f}")
