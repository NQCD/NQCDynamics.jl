

<p align="right">
  <a href="https://nqcd.github.io/NQCDynamics.jl/stable/">
    <img src="https://github.com/NQCD/NQCDLogo/blob/main/images/logo_with_text.png" alt="NQCDynamics.jl logo"
         title="NQCDynamics.jl" align="right" height="60"/>
  </a>
</p>

# NQCDynamics.jl

| **Documentation**                                     | **Build Status**                                |  **License**                     |
|:------------------------------------------------------|:----------------------------------------------- |:-------------------------------- |
| [![][docs-img]][docs-url] [![][ddocs-img]][ddocs-url] | [![][ci-img]][ci-url] [![][ccov-img]][ccov-url] | [![][license-img]][license-url]  |

[ddocs-img]: https://img.shields.io/badge/docs-dev-blue.svg
[ddocs-url]: https://nqcd.github.io/NQCDynamics.jl/dev/

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://nqcd.github.io/NQCDynamics.jl/stable/

[ci-img]: https://github.com/nqcd/NQCDynamics.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/nqcd/NQCDynamics.jl/actions/workflows/CI.yml

[ccov-img]: https://codecov.io/gh/NQCD/NQCDynamics.jl/branch/main/graph/badge.svg
[ccov-url]: https://codecov.io/gh/NQCD/NQCDynamics.jl

[license-img]: https://img.shields.io/github/license/NQCD/NQCDynamics.jl
[license-url]: https://github.com/NQCD/NQCDynamics.jl/blob/main/LICENSE

**Fast and flexible nonadiabatic molecular dynamics in Julia!**

-  🚗 **Fast:** uses [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) for efficient dynamics.
-  🪚 **Extensible:** plenty of room for more methods.
- ⚛️ **Transferable:** handles both simple models and atomistic systems.
- 👩‍🏫 **Helpful:** extended documentation with plenty of examples.

<p align="center">
<a href="https://nqcd.github.io/NQCDynamics.jl/stable/"><strong>Explore the NQCDynamics.jl docs 📚</strong></a>
</p>

---

With this package you can generate the initial conditions and perform the dynamics for your nonadiabatic dynamics simulations.
Tight integration with [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/)
makes the implementation of new methods relatively simple since we
build upon an already successful package providing a vast array of features.
We hope that the package will be of use to new students and experienced researchers alike, acting as a tool for learning and for developing new methods.

--- 

If you find this package to be useful, please cite the following paper:

J. Gardner et al., "NQCDynamics.jl: A Julia package for nonadiabatic quantum classical molecular dynamics in the condensed phase" <a href="https://pubs.aip.org/aip/jcp/article-abstract/156/17/174801/2841279/NQCDynamics-jl-A-Julia-package-for-nonadiabatic?redirectedFrom=fulltext"> J. Chem. Phys. 156, 174801 (2022) </a>
