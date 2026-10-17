# Hypersonic Cylinder

The bow shock in front of a cylinder at Mach 8.1, an Euler test of MOSE: HLLE, MUSCL with the van Albada limiter, RK3 with
a local time step, CFL 0.5, 5000 iterations. After 5000 iterations the two codes differ by 2.4e-2 on the temperature. The
reason is a property of the scheme on this case, present in both codes.

## The flow never reaches a steady state

The residual, the root mean square of the change of each variable per iteration, falls to its minimum at iteration 1551.
It then rises about six times, and stays there:

| Variant, run on MOSE | Density residual |
|---|---|
| as in the test: MUSCL, van Albada | 0.09–0.19 up to 50 000 iterations |
| MUSCL, minmod | 0.12 after 10 000 iterations |
| MUSCL, van Albada, CFL 0.25 | 0.097 after 10 000 iterations |
| first order | 1e-10 at iteration 6748, where the run stops on the residual threshold |

At first order the case converges to round-off. With the second-order reconstruction it does not, with either limiter
and at half the CFL number. The cause is the reconstruction at the captured shock, not the time step.

## Where and how

1 % of the cells holds 79 % of the squared residual. They are the cells of the bow shock in front of the stagnation line.
Probes in four of them show the density flipping at every iteration between about 1.3 kg/m³, the free-stream value, and
about 3 kg/m³. That is a period of two iterations, with a swing of 60–80 % of the mean. A probe inside the shock layer
swings by 3 %.

## What it means for the comparison

When the flow is deterministic, the two codes coincide. Once the flip-flop starts, each code differs from itself when the
CFL number changes by one part in 10^12, by 1–3e-2 after 5000 iterations: as much as from the other code. At first order
both codes converge, stop at the same iteration, and agree to 1.1e-13.
