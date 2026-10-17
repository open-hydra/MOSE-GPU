#!/usr/bin/env python3
"""Validate the chimera records of a MOSE bc file against the mesh.

MOSE builds a chimera ghost cell as the weighted sum of its donors' conserved
state, so a record is only meaningful if its weights sum to 1 and every donor
index exists. ATLAS has shipped files with neither property (donor lists that
sum to ~0, donors left null), and the solver's own guard is a runtime abort
deep into the run -- this checks the file up front.

    bc_chimera_check.py <bc file> <mesh.tec|.dat>

Reports, per (receiving block, face): record count, donor blocks, donors per
ghost layer, and the worst weight-sum error. Exits non-zero on any violation.
"""

import sys
from collections import defaultdict

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import tecmesh

CHIMERA_IDS = (102, 666)
ONE_PROPERTY = set([1, 4, 5, 6, 8, 9, 14, 22, 101, 103, 201, 410, 420, 501, 502, 1000, 1001, 999, 667]) \
    | set(range(301, 310)) | set(range(401, 409)) | {403, 407}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    bcfile, mesh = sys.argv[1], sys.argv[2]

    zones = tecmesh.zone_index(mesh)
    dims = {b + 1: tuple(n - 1 for n in z[1:4]) for b, z in enumerate(zones)}   # cell counts
    print("[mesh] %d blocks" % len(dims))

    stats = defaultdict(lambda: dict(n=0, donors=set(), lay=set(), werr=0.0, ndon=0))
    errors = []
    nrec = 0
    with open(bcfile) as fh:
        it = iter(fh)
        for line in it:
            p = line.split()
            if len(p) != 6:
                continue
            b, i, j, k, f, t = (int(x) for x in p)
            nrec += 1
            if t in CHIMERA_IDS:
                n1, n2 = (int(x) for x in next(it).split())
                groups = []
                for n in (n1, n2):
                    g = []
                    for _ in range(n):
                        d = next(it).split()
                        g.append((int(d[0]), int(d[1]), int(d[2]), int(d[3]), float(d[4])))
                    groups.append(g)
                s = stats[(b, f)]
                s["n"] += 1
                s["lay"].add((n1, n2))
                s["ndon"] += n1 + n2
                for gi, g in enumerate(groups, 1):
                    if not g:
                        errors.append("record %d (b=%d f=%d): ghost layer %d has no donors" % (nrec, b, f, gi))
                        continue
                    w = sum(d[4] for d in g)
                    s["werr"] = max(s["werr"], abs(w - 1.0))
                    if abs(w - 1.0) > 1e-9:
                        errors.append("record %d (b=%d,%d,%d,%d f=%d): layer %d weights sum to %.12f"
                                      % (nrec, b, i, j, k, f, gi, w))
                    for db, di, dj, dk, _ in g:
                        s["donors"].add(db)
                        if db not in dims:
                            errors.append("record %d: donor block %d does not exist" % (nrec, db))
                        else:
                            ni, nj, nk = dims[db]
                            if not (1 <= di <= ni and 1 <= dj <= nj and 1 <= dk <= nk):
                                errors.append("record %d: donor (%d,%d,%d,%d) outside block dims %s"
                                              % (nrec, db, di, dj, dk, dims[db]))
                        if db == b:
                            errors.append("record %d: block %d donates to itself" % (nrec, b))
                if len(errors) > 20:
                    break
            elif t in ONE_PROPERTY:
                next(it, None)

    if not stats:
        print("[chimera] no chimera records in this file")
        return 0

    print("[chimera] %-8s %-10s %-22s %-14s %s" % ("blk.face", "records", "donor blocks", "donors/layer", "max |sum w - 1|"))
    for (b, f) in sorted(stats):
        s = stats[(b, f)]
        print("          %-8s %-10d %-22s %-14s %.2e"
              % ("%d.%d" % (b, f), s["n"], sorted(s["donors"]),
                 sorted(s["lay"])[:3] if len(s["lay"]) <= 3 else "%d shapes" % len(s["lay"]),
                 s["werr"]))
    total = sum(s["n"] for s in stats.values())
    ndon = sum(s["ndon"] for s in stats.values())
    print("[chimera] %d records, %d donor entries (%.2f per record)" % (total, ndon, ndon / total))

    if errors:
        print("\n[FAIL] %d problem(s):" % len(errors))
        for e in errors[:20]:
            print("   ", e)
        return 1
    print("[OK] every donor exists and every ghost layer's weights sum to 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
