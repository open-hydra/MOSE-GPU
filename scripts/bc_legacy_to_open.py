#!/usr/bin/env python3
"""Convert a legacy MOSE bc.txt into the MOSE_open (ATLAS) 3-digit dialect.

MOSE_GPU reads both dialects (see src/lib/io/IO_BC_Open.f90); this converter exists
to build open-format twins of legacy cases, which is how the translation layer is
gated: legacy run and converted run must agree bit for bit.

Usage:
    bc_legacy_to_open.py <bc.txt> <out.txt> --nsc N [--nrans N]

The species/RANS counts are needed because the two dialects order the composition
columns differently: legacy writes (turb..., species...), MOSE_open writes
(species..., turb...).
"""

import argparse
import sys

# Direction columns: legacy alpha == beta == 0 means "use the face normal", which
# MOSE_open spells with the token below.
NORMAL_TOKEN = "normal"

# Legacy pamb sentinel meaning "inlet only; fall back to a wall on reversed flow".
PAMB_FORCE_INFLOW = -3.14

# Legacy mach-column encodings of the inflow modes.
MACH_MASSFLUX_T0 = -10.0
MACH_HRE_TANGENTIAL = -15.0


def fmt(values):
    """Format a property line the way ATLAS writes it (comma separated)."""
    out = []
    for v in values:
        out.append("%15s" % v if isinstance(v, str) else "%15.6E" % v)
    return ",".join(out) + ","


def is_number(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def direction(alpha, beta):
    """Legacy (alpha, beta) -> open direction tokens."""
    if alpha == 0.0 and beta == 0.0:
        return [NORMAL_TOKEN, NORMAL_TOKEN]
    return [alpha, beta]


def split_ci(ci, nsc, nrans):
    """Legacy composition order (turb, species) -> open order (species, turb)."""
    if len(ci) != nsc + nrans:
        raise ValueError("expected %d composition columns, got %d" % (nsc + nrans, len(ci)))
    return list(ci[nrans:]) + list(ci[:nrans])


def convert_inflow(props, nsc, nrans):
    """Legacy type 4 -> the matching MOSE_open inlet/ambient id."""
    mach = float(props[0])
    t0 = float(props[1])
    p0_tok = props[2]
    alpha, beta = float(props[3]), float(props[4])
    pamb, rel_fac = float(props[5]), float(props[6])
    ci = split_ci([float(x) for x in props[7:7 + nsc + nrans]], nsc, nrans)
    direc = direction(alpha, beta)

    if mach == MACH_HRE_TANGENTIAL:
        raise ValueError("legacy mach == -15 (HRE tangential injector) has no MOSE_open id")

    if mach < 0.0:
        g = float(p0_tok)
        if mach == MACH_MASSFLUX_T0:
            return 403, [t0, g] + direc + [rel_fac] + ci     # mass flux + total T
        if rel_fac != 1.0:
            raise ValueError("legacy mass-flux inlet with rel_fac != 1 is not a 404 (which "
                             "relaxes g towards the interior mass flux)")
        return 404, [t0, g] + direc + [rel_fac] + ci         # mass flux + static T

    if not is_number(p0_tok):                                 # time-varying p0 file
        return 402, [t0, p0_tok] + direc + [rel_fac] + ci
    p0 = float(p0_tok)

    if mach > 1.0:
        return 405, [mach, t0, p0] + direc + [rel_fac] + ci   # supersonic, static T/p
    if pamb == PAMB_FORCE_INFLOW:
        return 401, [t0, p0] + direc + [rel_fac] + ci         # stagnation inlet
    return 407, [t0, p0, pamb] + direc + [rel_fac] + ci       # ambient / far field


def convert(lines, nsc, nrans):
    out = []
    i = 0
    while i < len(lines):
        head = lines[i].split()
        if not head:
            i += 1
            continue
        b, ii, jj, kk, f, t = (int(x) for x in head[:6])
        i += 1

        def props():
            """Consume the property line and return its tokens."""
            nonlocal i
            toks = [x for x in lines[i].replace(",", " ").split() if x]
            i += 1
            return toks

        if t == 0:
            new_t, line = 0, None
        elif t == 1:
            new_t, line = 101, props()
        elif t == 2:
            new_t, line = 200, None
        elif t == 3:
            new_t, line = 300, None
        elif t == 11:
            new_t, line = 400, None
        elif t == 4:
            new_t, line = convert_inflow(props(), nsc, nrans)
        elif t == 5:
            new_t, line = 301, [float(props()[0]), 0.0, 0.0]
        elif t == 6:
            new_t, line = 302, [float(props()[0]), 0.0, 0.0]
        elif t == 8:
            p = props()
            new_t, line = 303, [float(p[0]), float(p[1]), 0.0]
        elif t in (9, 10, 12, 13):
            if t != 9:
                raise ValueError("pyrolysis material %d has no distinct MOSE_open id" % t)
            new_t, line = 304, [float(props()[0]), 0.0]
        elif t == 22:
            p = props()
            ci = split_ci([float(x) for x in p[7:7 + nsc + nrans]], nsc, nrans)
            # legacy nozzle: psub/psup live in the alpha/beta columns, in bar
            new_t = 420
            line = [float(p[1]), float(p[2]), float(p[3]) * 1e5, float(p[4]) * 1e5,
                    float(p[5]), float(p[6])] + ci
        elif t == 1001:
            new_t, line = 501, props()
        elif t == 1000:
            new_t, line = 103, props()
        elif t == 666:
            ni = props()
            n_donor = int(ni[0]) + int(ni[1])
            donors = [props() for _ in range(n_donor)]
            new_t, line = 102, None
            out.append("%8d%8d%8d%8d%8d%8d" % (b, ii, jj, kk, f, new_t))
            out.append("%8d%8d" % (int(ni[0]), int(ni[1])))
            for d in donors:
                out.append("%8d%8d%8d%8d%20s" % (int(d[0]), int(d[1]), int(d[2]), int(d[3]), d[4]))
            continue
        elif t == 667:
            p = lines[i].strip()
            i += 1
            new_t, line = 410, None
            parts = p.split()
            out.append("%8d%8d%8d%8d%8d%8d" % (b, ii, jj, kk, f, new_t))
            out.append(",".join(parts))
            continue
        elif t == 14:
            p = props()
            ci = [float(x) for x in p[5:5 + nsc]]
            # legacy: a, n, pRef, Taf, haf   ->   open: Taf, a, n, pRef, rhoGrain, SF
            new_t = 502
            line = [float(p[3]), float(p[0]), float(p[1]), float(p[2]), 0.0, 0.0] + ci
        else:
            raise ValueError("legacy BC type %d has no MOSE_open counterpart" % t)

        out.append("%8d%8d%8d%8d%8d%8d" % (b, ii, jj, kk, f, new_t))
        if line is not None:
            out.append(fmt(line) if not all(isinstance(x, str) for x in line) else " ".join(line))

    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--nsc", type=int, required=True, help="number of chemical species")
    ap.add_argument("--nrans", type=int, default=0, help="number of RANS variables")
    args = ap.parse_args()

    with open(args.infile) as fh:
        lines = fh.read().splitlines()

    try:
        out = convert(lines, args.nsc, args.nrans)
    except ValueError as err:
        sys.exit("[ERROR] %s" % err)

    with open(args.outfile, "w") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
