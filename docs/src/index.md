# Introduction

Welcome to the documentation for NonadiabaticMolecularDynamics, 
a package for performing nonadiabatic molecular dynamics simulations.
The documentation covers both how to use the existing code and describes the
intricacies of the implementations, hoping to make further contributions as simple as possible.

The goal of the package is to provide highly performant code that is easily extensible and
transferrable between model systems and realistic atomic systems.
The dynamics methods are built on top of the `DifferentialEquations.jl` ecosystem which
provides an intuitive interface for implementing different nonadiabatic dynamics methods.

The existing methods are :

- Classical molecular dynamics (MD)
- Classical Langevin dynamics (LD)
- Fewest-switches surface hopping (FSSH)
- Molecular dynamics with electronic friction (MDEF)
- Ring polymer molecular dynamics (RPMD)
- Nonadiabatic ring polymer molecular dynamics (NRPMD)

Given that mixed quantum-classical and semiclassical methods usually require averaging
over many trajectories, we provide a few options for generating initial conditions for
the simulations:

- Metropolis-Hastings Markov chain Monte Carlo thermal sampling
- Langevin dynamics

## Installation

The easiest way to install the package is to the install the `NonadiabaticRegistry`:
```julia
pkg> registry add "https://github.com/maurergroup/NonadiabaticRegistry"
```

Then you can easily add the package with:
```julia
pkg> add NonadiabaticMolecularDynamics
```

Included in the registry alongside the main package are a few others that provide extra
models and add-ons. Here is an overview of the currently existing packages included in
the registry:

```@raw html
<img src="./assets/registry.png">
```
