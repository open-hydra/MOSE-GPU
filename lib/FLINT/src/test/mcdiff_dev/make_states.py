#!/usr/bin/env python3
"""gpu-port-mcdiff unit test: states for test_co_DS_dev.
Real states: THORBL2 pocket slices (x 12 and x 25 mm, 1800 x 68 cells, MR fields +396/576/612/648/684 us), rho_i, T, p.
Edge cases appended: pure H2 / N2 / O2 / H (near-pure fallback), air with exact zeros, T below 1 K and above the table
(clamp), low and high pressure, and random compositions (Dirichlet) over 200-4500 K and 1e4-3e6 Pa.
usage: make_states.py out.bin pocket_mr*.npz"""
import sys, numpy as np
out, files = sys.argv[1], sys.argv[2:]
R, T, P = [], [], []
for f in files:
    d = np.load(f)
    for st in ("x12", "x25"):
        R.append(np.array([d[f"{st}_rho({s})"].ravel() for s in range(1, 8)], dtype=np.float64))
        T.append(d[f"{st}_T"].ravel().astype(np.float64)); P.append(d[f"{st}_p"].ravel().astype(np.float64))
R, T, P = np.concatenate(R, 1), np.concatenate(T), np.concatenate(P)
nreal = T.size
print(f"real states {nreal}: T {T.min():.1f}-{T.max():.1f} K, p {P.min():.0f}-{P.max():.0f} Pa, "
      f"rho_i min {R.min():.3e}, zeros {(R == 0).sum()}, negatives {(R < 0).sum()}")
E_R, E_T, E_P = [], [], []
def add(r, t, p): E_R.append(np.array(r, float)); E_T.append(float(t)); E_P.append(float(p))
for k in (0, 2, 3, 6):                                    # pure O2, H2, H, N2
    r = np.zeros(7); r[k] = 0.5; add(r, 1200.0, 1.3e5)
air = np.array([0.232, 0, 0, 0, 0, 0, 0.768]) * 1.2
for t in (0.3, 1.0, 1.7, 299.99, 4999.5, 5000.0, 5000.4, 7000.0, 12000.0):   # clamp below/above the table
    add(air, t, 1.0e5)
for p in (1.0e3, 1.0e4, 5.0e6):
    add(air, 800.0, p)
rng = np.random.default_rng(20260925)
Y = rng.dirichlet(np.full(7, 0.3), 10000)
for y, t, p in zip(Y, rng.uniform(200, 4500, 10000), 10 ** rng.uniform(4, np.log10(3e6), 10000)):
    add(y * 0.5, t, p)
R = np.concatenate([R, np.array(E_R).T], 1); T = np.concatenate([T, E_T]); P = np.concatenate([P, E_P])
n = T.size
with open(out, "wb") as f:
    np.array([n], np.int32).tofile(f)
    np.asfortranarray(R).T.ravel(order="C").astype("<f8").tofile(f)   # rhoi(7,n) column-major = R.T row-major
    T.astype("<f8").tofile(f); P.astype("<f8").tofile(f)
print(f"written {out}: {n} states ({nreal} real + {n - nreal} edge/random)")
