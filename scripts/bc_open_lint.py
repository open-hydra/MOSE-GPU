#!/usr/bin/env python3
"""Walk MOSE_open (ATLAS) bc files with the same record structure MOSE_GPU's
reader uses, and report what they contain.

Purpose: check a real ATLAS case against the dual-mode reader without running
the solver. A file that does not walk cleanly here would desynchronise the
Fortran prescan too (the boundary-face count would not match the IC).

Usage:  bc_open_lint.py <bc file> [more files ...]
"""

import sys

# Records with no property line at all (agrees with MOSE_open's Check_BC and
# with Skip_BC_Record in src/lib/io/IO_BC_Open.f90).
NO_PROPERTY = {0, 100, 200, 300, 400, 409}
# Records with one property line.
ONE_PROPERTY = set([101, 103, 201, 410, 420, 501, 502]) | set(range(301, 310)) | set(range(401, 409))
# 102 (chimera) is variable length: an "ni" line plus one line per donor.

# Types MOSE_GPU refuses or restricts (kept in sync with IO_BC_Open.f90).
RESTRICTED = {
    404: "exact only with rel_fac = 1",
    405: "requires Mach > 1; 'normal' direction only on x-aligned faces",
    409: "not supported",
}


def walk(path):
    lines = open(path, errors="replace").read().splitlines()
    types, i, nrec = {}, 0, 0
    props = {}
    while i < len(lines):
        head = lines[i].split()
        if not head:
            i += 1
            continue
        if len(head) < 6:
            raise ValueError("line %d is not a record header: %r" % (i + 1, lines[i]))
        t = int(head[5])
        i += 1
        nrec += 1
        types[t] = types.get(t, 0) + 1
        if t == 102:
            n = lines[i].split()
            i += 1 + int(n[0]) + int(n[1])
        elif t in ONE_PROPERTY:
            props.setdefault(t, set()).add(len([x for x in lines[i].replace(",", " ").split() if x]))
            i += 1
        elif t not in NO_PROPERTY:
            raise ValueError("unknown BC type %d at line %d" % (t, i))
    return nrec, types, props


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    rc = 0
    for path in sys.argv[1:]:
        try:
            nrec, types, props = walk(path)
        except (ValueError, IndexError) as err:
            print("%s: [ERROR] %s" % (path, err))
            rc = 1
            continue
        dialect = "MOSE_open" if any(t >= 100 for t in types) else "legacy"
        print("%s: %s, %d records" % (path, dialect, nrec))
        for t in sorted(types):
            note = ""
            if t in props:
                note = " (%s property columns)" % "/".join(str(c) for c in sorted(props[t]))
            if t in RESTRICTED:
                note += "  [MOSE_GPU: %s]" % RESTRICTED[t]
            print("    %4d : %8d%s" % (t, types[t], note))
    sys.exit(rc)


if __name__ == "__main__":
    main()
