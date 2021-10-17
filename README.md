

<p align="right">
  <a href="https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/">
    <img src="https://github.com/NQCD/NQCDLogo/blob/main/logo.png" alt="NonadiabaticMolecularDynamics.jl logo"
         title="NonadiabaticMolecularDynamics.jl" align="right" height="60"/>
  </a>
</p>

# NonadiabaticMolecularDynamics.jl

[![CI](https://github.com/nqcd/NonadiabaticMolecularDynamics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/nqcd/NonadiabaticMolecularDynamics.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/NQCD/NonadiabaticMolecularDynamics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NQCD/NonadiabaticMolecularDynamics.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/)

**Fast and flexible nonadiabatic molecular dynamics in Julia!**

-  🚗 **Fast:** uses [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) for efficient dynamics.
-  🪚 **Extensible:** plenty of room for more methods.
- ⚛️ **Transferable:** handles both simple models and atomistic systems.
- 👩‍🏫 **Helpful:** extended documentation with plenty of examples.

<p align="center">
<a href="https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/"><strong>Explore the NonadiabaticMolecularDynamics.jl docs 📚</strong></a>
</p>

---

With this package you can generate the initial conditions and perform the dynamics for your nonadiabatic dynamics simulations.
Tight integration with [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)
makes the implementation of new methods relatively simple since we
build upon an already successful package providing a vast array of features.
We hope that the package will be of use to new students and experienced researchers alike, acting as a tool for learning and for developing new methods.
