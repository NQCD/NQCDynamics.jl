# [Introduction](@id dynamicssimulations)

Performing dynamics simulations is at the core of this package's functionality
(as you might have guessed from the name).
This section of the documentation will describe how to perform dynamics simulations,
building on the introduction from [Getting started](@ref).

Since we use [DifferentialEquations](https://diffeq.sciml.ai/stable/)
to perform the dynamics, it is most natural
to split up the system parameters from the dynamics variables.
This manifests itself as two separate data types: the [`Simulation`](@ref), and the
[`DynamicsVariables`](@ref).

!!! info

    If you intend to implement a new dynamics method, we recommend reading
    [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) to understand
    more deeply how this package works.

The [`Simulation`](@ref) holds all the static information about the system: the atoms,
the model, the temperature, the cell and the dynamics method.

```@example dynamics
using NQCDynamics # hide
atoms = Atoms(2000) # Single atom with mass = 2000 a.u.
sim = Simulation{Ehrenfest}(atoms, TullyModelOne(); temperature=0, cell=InfiniteCell())
```
Here we have initialised the simulation parameters, including the default temperature and cell explicitly.
`sim` takes the place of the `p` parameter seen throughout [DifferentialEquations](https://diffeq.sciml.ai/stable/).

For [DifferentialEquations](https://diffeq.sciml.ai/stable/) to allow for a wide variety of solvers, 
the input arrays [`DynamicsVariables`](@ref) must be [`AbstractArray`](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)s.
In nonadiabatic dynamics simulations, we usually have different groups of variables that behave in particular ways.
For example: for mapping variable methods we have positions, velocities, and two sets of mapping variables representing
the electronic degrees of freedom.

For this purpose we use the [`ComponentVector`](https://github.com/jonniedie/ComponentArrays.jl), which allows us
to arbitrarily partition the variables into their subgroups.
This allows us to keep all the variables in a single array as required by
[DifferentialEquations](https://diffeq.sciml.ai/stable/),
whilst still having them partitioned for convenient computation and readable code.

```@example dynamics
v0 = hcat(10) / 2000
r0 = hcat(-5)
u0 = DynamicsVariables(sim, v0, r0, PureState(1))
```

Since each dynamics method has a different set of variables, each method implements
[`DynamicsVariables(sim, args...)`](@ref DynamicsVariables), which will convert the input into the correct
structure.
This helps to ensure each method follows a similar workflow, making it easy to switch between different methods.
The output of this function takes the place of the `u` argument seen throughout
[DifferentialEquations](https://diffeq.sciml.ai/stable/).

With both the [`Simulation`](@ref) and [`DynamicsVariables`](@ref) in hand,
the central function is [`run_dynamics`](@ref) which allows us to perform a single dynamics trajectory.
[`run_dynamics`](@ref) takes the simulation parameters `sim` and the initial conditions `u0`, along with a time span `tspan`
that the trajectory will cover.

```@example dynamics
tspan = (0.0, 2000.0)
run_dynamics(sim, tspan, u0; output=OutputDynamicsVariables, dt=1.0)
```

The output is a dictionary containing entries for `:Time` and our requested output quantity. 
Output is a required keyword and the code will error unless at least one quantity is specified.
By passing a `Tuple` to the `output` keyword argument we can ask for multiple quantities.

```@example dynamics
out = run_dynamics(sim, tspan, u0; output=(OutputPosition, OutputAdiabaticPopulation))
```
The quantities that are available are listed [here](@ref NQCDynamics).
More quantities can be added by defining new functions with the signature `f(sol, i)`.
The first argument is the `DifferentialEquations.jl` solution object and the second is the trajectory index.

This time we can see that the output contains the two quantities that we asked for.

```@example dynamics
using Plots
plot(out, :OutputPosition)
```

```@example dynamics
plot(out, :OutputAdiabaticPopulation)
```

!!! note

    Here we have used a special plot recipe that will automatically plot any quantity against time.
    This is useful when investigating the results of a single trajectory.

All of the dynamics methods work in a similar way. For details on a specific method along with examples,
please see the method specific page in the sidebar under `Dynamics methods`.
