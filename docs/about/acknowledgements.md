# Acknowledgements

MOSE-GPU is built upon several projects.

## MOSE

**Multipurpose Open Solver for ideal-gas Equations**

- **Repository:** [github.com/open-hydra/MOSE](https://github.com/open-hydra/MOSE)
- **License:** GPL v3.0

MOSE-GPU is the GPU port of MOSE: its equations, numerical methods, input files and test cases are those of MOSE.

## FLINT

**Fortran Library for INtegrated Thermochemistry**

- **Repository:** [github.com/MarcoGrossi92/FLINT](https://github.com/MarcoGrossi92/FLINT)
- **License:** GPL v3.0

Thermodynamics, transport and chemistry. MOSE-GPU includes a version of FLINT with the device integrators and the device
thermodynamic and transport routines.

## OSlo

- **Repository:** [github.com/MarcoGrossi92/OSlo](https://github.com/MarcoGrossi92/OSlo)
- **License:** GPL v3.0

The ODE solver library used by FLINT.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

Reading and writing of Tecplot, VTK and Plot3D files, including the Tecplot binary format through TecIO.

## FiNeR

**Fortran INI ParseR and generator**

- **Repository:** [github.com/szaghi/FiNeR](https://github.com/szaghi/FiNeR)
- **License:** GPL v3.0

## Documentation Tools

### MkDocs and Material for MkDocs

- **Websites:** [mkdocs.org](https://www.mkdocs.org), [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
- **Contribution:** documentation site generation and theme
