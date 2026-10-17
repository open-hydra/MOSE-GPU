# 🔥 FLINT – Fortran Library for INTegrated Thermochemistry

FLINT is a lightweight, modular Fortran library for accessing thermodynamic, transport, and chemical source term data. It optionally uses Cantera as a backend, making it ideal for integration with CFD solvers and reactive flow applications.

## 🧱 Why FLINT?

FLINT provides a clean interface layer that lets you:

- Focus on pure thermochemistry, without getting bogged down by unrelated CFD complexities
- Extend or override functionality without rewriting the entire chemistry module
- Build and test your chemical models in a modular, maintainable way
- Use Cantera to compute the chemical source terms when needed

Like the stone it’s named after, FLINT is small but powerful — it provides the spark that ignites complex thermochemical simulations.

## ✨ Features

- **Supports species thermodynamics, chemical source terms, and transport coefficients**
- **Designed for coupling with CFD solvers**, especially for combustion and reacting flows
- **Lightweight and fast**, with minimal dependencies
- **Modular architecture** - plug in your own chemistry or transport data sources
- **Optional Cantera integration** for advanced kinetic models
- **Built-in library of chemical models** tested and verified
- **Source code generator** for converting a yaml file to an F90 source code to be appended to the FLINT library
- **General routine** for specific kinetic models not contained into FLINT library 

## 📦 Getting Started

```fortran
use FLINT_Load_ThermoTransport
use FLINT_Load_Chemistry

err = read_idealgas_thermo(folder='WD/INPUT')
err = read_idealgas_transport(folder='WD/INPUT')
err = read_chemistry(folder='WD/INPUT', mech_name=mech_name)
```

## 🛠️ Contributing: Adding a New Mechanism

Follow these steps to integrate a new chemical mechanism into FLINT:

0. Build FLINT againts Sundials and Cantera. If it is not possible, bypass steps 5 and 6 and use **KAnt** to generate a benchmark.

1. ✅ Verify the YAML Mechanism
- Ensure your `<mech>.yaml` file is valid and behaves as expected.
- Use external tools like **KAnT** or **Cantera** for validation.

2. 📁 Create a Dedicated Test Folder
- In the main test directory, create a subfolder named `<mech>`.
- Build and store all relevant input, output, and test artifacts in this folder.

3. ⚙️ Generate Source Code Using YTF
- Use `utils/YTF` (YAML-to-Fortran) to convert the `<mech>.yaml` file into Fortran-compatible source code.
- Modify or adapt the generated code as needed to integrate with FLINT.

4. 🧪 Update Fortran Batch Tests
- Add a new test section for `<mech>` in the Fortran test suite.
- Follow the structure of existing mechanism tests.

5. 🧪 Update C++ Batch Tests
- Add a corresponding test case for `<mech>` in the C++ batch tests.
- Follow the structure of existing mechanism tests.

6. 🔍 Verify Consistency
- Run both Fortran and C++ tests.
- Confirm that results are consistent across both implementations.

---

> More detailed setup and examples coming soon.
