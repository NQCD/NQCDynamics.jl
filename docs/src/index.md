# Introduction

Welcome to the documentation for NonadiabaticMolecularDynamics, 
a package for performing nonadiabatic molecular dynamics simulations.
The documentation covers both how to use the existing code and describes the
intricacies of the implementations, hoping to make further contributions as simple as possible.

### Objectives

- Achieve high performance along with good readability, extensibility, maintainability
- Handle both simple models and high-dimensional systems
- Highlight the advantages of Julia in the field of nonadiabatic dynamics
- Encourage code sharing and reuse within the nonadiabatic dynamics community

Reproducibility is a pressing issue in the field of theoretical chemistry and physics as often studies either do not attempt to provide all necessary data or code for full reproducibility of the work. This can lead to difficulties when attempting to better understand the theory and implementation of the method and makes it difficult for students not only to learn existing models and theories, but also to improve and extend these. 
This project provides implementations for existing dynamics methods along with
a framework that can be used for future research with the goal of encouraging greater
code sharing and reuse within the nonadiabatic dynamics community.

### Features

Here we provide a list of currently implemented features of the code.
We encourage contributions and implementations of methods. To do so, please open up an
issue/pull request on Github!

#### Dynamics methods

- [Classical molecular dynamics](@ref classical-dynamics)
- [Classical Langevin dynamics](@ref langevin-dynamics)
- [Fewest-switches surface hopping (FSSH)](@ref fssh-dynamics)
- [Molecular dynamics with electronic friction (MDEF)](@ref mdef-dynamics)
- [Ring polymer molecular dynamics (RPMD)](@ref rpmd-dynamics)
- [Nonadiabatic ring polymer molecular dynamics (NRPMD)](@ref nrpmd-dynamics)
- [Ring polymer surface hopping (RPSH)](@ref rpsh-dynamics)
- [Ehrenfest molecular dynamics](@ref ehrenfest-dynamics)

#### Generating initial conditions

- [Thermal Metropolis-Hastings Monte Carlo](@ref mhmc-sampling)
- [Thermal Hamiltonian Monte Carlo](@ref hmc-sampling)
- [Langevin dynamics](@ref langevin-sampling)
- [Semiclassical EBK quantisation](@ref ebk-sampling)

### Dynamics with `DifferentialEquations.jl`

The [`DifferentialEquations`](https://diffeq.sciml.ai/stable/) ecosystem from the
[SciML organisation](https://github.com/SciML/) provides a large library of integration
algorithms along with a simple interface for implementing new algorithms that can be tailored
for specific nonadiabatic dynamics methods.
Further, they provide helpful utilities for introducing discontinuities through the 
[callback interface](https://diffeq.sciml.ai/stable/features/callback_functions/#Using-Callbacks)
or handling many trajectories at once to obtain ensemble averaged observables with
the [ensemble interface](https://diffeq.sciml.ai/stable/features/ensemble/).
We can take advantage of these utilities by basing our dynamics setup on this framework
which significantly simplifies the implementation of new methods.

### Installation

#### 1. Install Julia
Download and install the current stable release from the [Julia website](https://julialang.org/downloads/).
For most platforms `julia` is provided as a precompiled binary and do not require any installation procedure. However, you need to specify the path to julia or create a symbolic link to the executable that is in your systempath. 

#### 2. Install the `NQCDRegistry`
Since the package is not included in the default registry (`General`), we must first
install the `NQCDRegistry`.
This gives access to the core `NonadiabaticMolecularDynamics` package along with some dependencies
and add-ons.
First, enter the Julia REPL by executing `julia` from the command line.
Then press `]` to enter `pkg` mode. The prompt should change from `julia>` to `pkg>`.
Install the registry directly from Github with: 
```julia-repl
pkg> registry add "https://github.com/NQCD/NQCDRegistry"
```

#### 3. Install the package
Now that the registry has been added, the package can be installed in the same way as any other registered Julia package:
```julia-repl
pkg> add NonadiabaticMolecularDynamics
```

#### 4. Use the package!
```julia-repl
julia> using NonadiabaticMolecularDynamics
```

To check the package has been installed correctly and everything is working, you can execute the tests
with:
```julia-repl
pkg> test NonadiabaticMolecularDynamics
```
Alternatively, you can proceed directly the next section for a walkthrough of some basic functionality.

### How to use this documentation

The first page to read is the [Getting started](@ref) section which walks through all the ingredients
needed to perform a conventional classical molecular dynamics simulation.
After this, the reader is free to explore at their leisure since everything else builds directly
upon sections from the [Getting started](@ref) page.

### Package ecosystem

Included in the `NQCDRegistry` alongside the main package are a few others that provide extra
models and add-ons. Here is an overview of the currently existing packages included in
the registry:

```@raw html
<img src="./assets/registry.png">
```
