# Introduction

Performing dynamics simulations is at the core of this package's functionality
(as you might have guessed from the name).
This section of the documentation will describe how to perform dynamics simulations,
building on the simple introduction from [`Getting started`].

Since we use [DifferentialEquations](https://diffeq.sciml.ai/stable/)
to perform the dynamics, it is most natural
to split up the system parameters from the dynamics variables.
This manifests itself as two separate data types: the [`Simulation`](@ref), and the
[`DynamicsVariables`](@ref).

!!! info

    An basic understanding of [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) helps to understand
    how this package works, especially if you intend to implement a new dynamics method.

The [`Simulation`](@ref) holds all the static information about the system: the atoms,
the model, the temperature, the cell and the dynamics method.
This takes the place of the `p` parameter seen throughout [DifferentialEquations](https://diffeq.sciml.ai/stable/).

```@example dynamics
using NonadiabaticMolecularDynamics # hide
atoms = Atoms(2000) # Single atom with mass = 2000 a.u.
sim = Simulation{Ehrenfest}(atoms, TullyModelOne(); temperature=0, cell=InfiniteCell())
```
Here we have initialised the simulation parameters, including the default temperature and cell explicitly.

The [`DynamicsVariables`](@ref) must be [`AbstractArray`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)s,
which means behave similarly to a standard array of numbers.
This is a requirement of [DifferentialEquations](https://diffeq.sciml.ai/stable/) which allows for the wide
variety of solvers to be used. 
In nonadiabatic dynamics simulations, we usually have different groups of variables that behave in particular ways.
For example: for mapping variable methods we have positions, velocities, and two sets of mapping variables representing
the electronic degrees of freedom.

For this purpose we use the [`ComponentVector`](https://github.com/jonniedie/ComponentArrays.jl) which allows us
to arbitrarily partition the variables into their subgroups.
This allows us to keep all the variables in a single array as required by
[DifferentialEquations](https://diffeq.sciml.ai/stable/),
whilst still having them partitioned for convenient computation and readable code.

```@example dynamics
v0 = hcat(10) / 2000
r0 = hcat(-5)
u0 = DynamicsVariables(sim, v0, r0, 1)
```

Since each dynamics method has a different set of variables, each method implements
[`DynamicsVariables(sim, args...)`](@ref DynamicsVariables) which will convert the input into the correct
structure.
This helps to ensure each method follows a similar workflow making for easy switching between different methods.
The output of this function takes the place of the `u` argument seen throughout
[DifferentialEquations](https://diffeq.sciml.ai/stable/).

With both the [`Simulation`](@ref) and [`DynamicsVariables`](@ref) in hand,
the central function is [`run_trajectory`](@ref) which allows us to perform a single dynamics trajectory.

```@example dynamics
run_trajectory(u0, (0.0, 2000.0), sim)
```

By default, we can see that the output contains both `t` and `u`. These are the time and dynamics variables
respectively.
By passing a `Tuple` to the `output` keyword argument we can ask for specific quantities.

```@example dynamics
out = run_trajectory(u0, (0.0, 2000.0), sim; output=(:position, :adiabatic_population))
```

This time we can see that the output contains only the quantities that we asked for.

```@example dynamics
using Plots
plot(out, :position)
```

```@example dynamics
plot(out, :adiabatic_population)
```

!!! note

    Here we have used a special plot recipe that will automatically plot any quantity against time.
    This is useful when investigating the results of a single trajectory.

All of the dynamics methods work in a similar way. For details on a specific method along with examples,
please see the method specific page in the sidebar under `Dynamics methods`.
