```@setup logging
@info "Expanding src/index.md..."
start_time = time()
```
# Introduction

Welcome to the documentation for NQCDynamics, 
a package for performing nonadiabatic molecular dynamics simulations.
The documentation covers both how to use the existing code and describes the
intricacies of the implementations, hoping to make further contributions as simple as possible.

### Objectives

- Achieve high performance along with good readability, extensibility, maintainability
- Handle both simple models and high-dimensional systems
- Highlight the advantages of Julia in the field of nonadiabatic dynamics
- Encourage code sharing and reuse within the nonadiabatic dynamics community

Reproducibility is a pressing issue in the fields of theoretical chemistry and computational materials physics as many studies do not provide all necessary data or code for full reproducibility of the work. This can lead to difficulties when attempting to better understand the theory and implementation of the method and makes it difficult for students to learn existing models and theories, but also to improve upon them. 
This project provides implementations for existing dynamics methods along with a framework that can be used for future research with the goal of encouraging greater code sharing and reuse within the nonadiabatic dynamics community.

### Features

Here we provide a list of some of the currently implemented features of the code. 
We encourage contributions and implementations of methods. To do so, please open up an
issue/pull request on Github!

#### Dynamics methods

- [Classical molecular dynamics](@ref classical-dynamics)
- [Classical Langevin dynamics](@ref langevin-dynamics)
- [Ehrenfest molecular dynamics](@ref ehrenfest-dynamics)
- [Fewest-switches surface hopping (FSSH)](@ref fssh-dynamics)
- [Independent Electron Surface Hopping (IESH)](@ref iesh-dynamics)
- [Molecular dynamics with electronic friction (MDEF)](@ref mdef-dynamics)
- [Ring polymer molecular dynamics (RPMD)](@ref rpmd-dynamics)
- [Nonadiabatic ring polymer molecular dynamics (NRPMD)](@ref nrpmd-dynamics)
- [Ring polymer surface hopping (RPSH)](@ref rpsh-dynamics)


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
For most platforms `julia` is provided as a precompiled binary and does not require any installation procedure. However, you need to specify the path to julia or create a symbolic link to the executable that is in your systempath. 

#### 2. OPTIONAL - Install the `NQCRegistry`
The `NQCDynamics.jl` package is included in the default registry (`General`) and so can be installed according to step 3. 
Most ussers can simply skip this step. For development purposes, it might be beneficial to add the NQCRegistry.

The user may also install `NQCDynamics.jl` and its dependencies from the `NQCRegistry` which requires manually adding.
First, enter the Julia REPL by executing `julia` from the command line.
Then press `]` to enter `pkg` mode. The prompt should change from `julia>` to `pkg>`.
Install the registry directly from Github with: 
```julia-repl
pkg> registry add "https://github.com/NQCD/NQCRegistry"
```

!!! warning

    If this is the first time you're using Julia there's a chance that the
    [General registry](https://github.com/JuliaRegistries/General) will not have been
    installed. Run `pkg> registry status` to view the installed registries.
    If `General` is not present, run `pkg> registry add General` before proceeding
    to the next step.

#### 3. Install the package
The package can be installed in the same way as any other registered Julia package:
```julia-repl
pkg> add NQCDynamics
```

#### 4. Use the package!
```julia-repl
julia> using NQCDynamics
```
You are now free to proceed to the next section and learn how to use the package.
If you would like you can complete step 5 to double check your installation.

#### 5. Run the tests (optional)

To check the package has been installed correctly and everything is working,
you can execute the tests with:
```julia-repl
pkg> test NQCDynamics
```

!!! warning

    The tests use some extra functionality from the
    [JuliaMolSim registry](https://github.com/JuliaMolSim/MolSim)
    which can be added directly from Github with
    `pkg> registry add "https://github.com/JuliaMolSim/MolSim"`.
    Without this, the tests will not run successfully.

### How to use this documentation

The first page to read is the [Getting started](@ref) section which walks through all the ingredients
needed to perform a conventional classical molecular dynamics simulation.
After this, the reader is free to explore at their leisure since everything else builds directly
upon sections from the [Getting started](@ref) page.

### How to cite
If you find this package to be useful, please cite the following paper:

- J. Gardner et al., ["NQCDynamics.jl: A Julia package for nonadiabatic quantum classical molecular dynamics in the condensed phase", J. Chem. Phys. 156, 174801 (2022)](https://pubs.aip.org/aip/jcp/article-abstract/156/17/174801/2841279/NQCDynamics-jl-A-Julia-package-for-nonadiabatic?redirectedFrom=fulltext)

When using the [independent-electron surface hopping (IESH)](@ref iesh-dynamics) implementation, please cite:

- J. Gardner et al., ["Efficient implementation and performance analysis of the independent electron surface hopping method for dynamics at metal surfaces", J. Chem. Phys. 158, 064101 (2023)](https://pubs.aip.org/aip/jcp/article/158/6/064101/2874963/Efficient-implementation-and-performance-analysis)

Benchmarking information of multiple mixed-quantum-classical (MQC) dynamics methods implemented within NQCD is given in the following paper:

- J. Gardner et al., ["Assessing Mixed Quantum-Classical Molecular Dynamics Methods for Nonadiabatic Dynamics of Molecules on Metal Surfaces", J. Phys. Chem. C 127, 15257 (2023)](https://pubs.acs.org/doi/10.1021/acs.jpcc.3c03591)

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
