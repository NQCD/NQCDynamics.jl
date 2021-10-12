# NonadiabaticMolecularDynamics.jl

[![CI](https://github.com/nqcd/NonadiabaticMolecularDynamics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/nqcd/NonadiabaticMolecularDynamics.jl/actions/workflows/CI.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/dev/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/stable/)

This is a package for performing nonadiabatic molecular dynamics simulations for both simple models and full-dimensional systems using the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) framework.

## Included functionality

This package provides ways to generate initial conditions for the dynamics simulations and perform the simulations
using DifferentialEquations.jl.
The tight integration with DifferentialEquations.jl makes the implementation of new methods relatively simple,
making it necessary only to define the time-derivatives of the dynamical variables.
This makes prototyping new methods as simple as possible.

Further, the model interface has been written in a dimension-agnostic way such that the same code can be used
for simple models and for full dimensional systems.

### Sampling initial conditions
- Monte carlo sampling
- Thermal Langevin dynamics
- Nonequilibrium quantisation of diatomic molecules

### Dynamics methods
- Fewest-switches surface hopping
- Ehrenfest molecular dynamics
- Nonadiabatic ring polymer molecular dynamics
- Ring polymer surface hopping
- Ehrenfest ring polymer molecular dynamics
- Molecular dynamics with electronic friction

For further details please consult the [documentation](https://nqcd.github.io/NonadiabaticMolecularDynamics.jl/dev/).
