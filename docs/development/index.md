# Developer Guide

| Page | Content |
|---|---|
| [Code Structure](structure.md) | where things live, and the two forms of the state |
| [Testing](testing.md) | the regression table and the comparison with MOSE |
| [Build Flags](build-flags.md) | the reference build and its compile-time options |
| [Contributing](contributing.md) | how to propose a change |

Two rules follow from the purpose of MOSE-GPU, which is to give the results of MOSE on GPUs:

1. **Physics and numerics follow MOSE.** A change to a method goes into MOSE first, or into both codes together, with
   the comparison of [Testing](testing.md).
2. **One code path.** There is no CPU fallback and no alternative kernel behind a switch. A faster kernel replaces the old
   one when the regression table is unchanged.
