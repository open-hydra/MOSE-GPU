#!/usr/bin/env python3
"""gpu-port-mcdiff unit test inputs: thermo.dat of the THOR runs cut at T <= TMAX (the 7-species diffusion.dat stops at 5000 K and
read_idealgas_diffusion requires the same Tmax as the thermo table). usage: make_input.py thermo_in thermo_out TMAX"""
import sys, re
src, dst, tmax = sys.argv[1], sys.argv[2], int(sys.argv[3])
out, keep = [], 0
for line in open(src):
    if line.startswith("ZONE"):
        keep = 0; out.append(line); continue
    m = re.match(r"\s*I=\s*(\d+)", line)
    if m:
        out.append(re.sub(r"I=\s*\d+", f"I={tmax}", line)); keep = tmax; continue
    if keep > 0 and re.match(r"\s*[0-9.]", line):
        out.append(line); keep -= 1; continue
    if not re.match(r"\s*[0-9.]", line): out.append(line)
open(dst, "w").writelines(out)
