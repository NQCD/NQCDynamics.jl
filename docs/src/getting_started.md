```@setup logging
@info "Expanding src/getting_started.md..."
start_time = time()
```
# Getting started

To get started with the package we can identify the necessary ingredients to
perform a simple classical dynamics simulation and walk through how
to set up the simulation.

```@setup started
using NQCDynamics
```

### Atoms

First, we must define the particles in the simulation.
For this purpose we provide the [`Atoms`](@ref NQCBase.Atoms) 
type which will contain the symbols, atomic numbers and masses for our atoms.
Technically these need not be actual atoms and be a generic particle.

If using real atoms, then they can be constructed using the chemical symbols
as a `Vector` of Julia's `Symbol` types, a `Vector{Symbol}`:
```@repl started
Atoms([:H, :C])
```
You can see that this contains two atoms labelled by their atomic numbers
with their masses in atomic units.

!!! note "Atomic units"

    Internally [atomic units](https://en.wikipedia.org/wiki/Hartree_atomic_units) are used
    for all quantities. This makes things simple when performing nonadiabatic dynamics.
    [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) and
    [UnitfulAtomic.jl](https://github.com/sostock/UnitfulAtomic.jl) can be used to help
    with unit transformations, and many functions will directly accept
    Unitful quantities and handle the conversions for you.

Alternatively, if not using real atoms, [`Atoms`](@ref NQCBase.Atoms)
can be created using a `Vector{<:Real}` where the provided numbers are the masses of the
particles.
```@repl started
Atoms([1, 2, 3, 4, 5, 6])
```

A more detailed look into the [`Atoms`](@ref NQCBase.Atoms) type along
with a description of how to save and load structures can be found
[here](@ref atoms).

### Representing atomic positions and velocities 

This package chooses to separate the dynamical variables from the static atomic parameters
included in the [`Atoms`](@ref NQCBase.Atoms) type.
This allows us to easily interface with other numerical packages like
[DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) and
[AdvancedMH.jl](https://github.com/TuringLang/AdvancedMH.jl).
As such, both positions and velocities are represented using Julia's standard `Array`
type, specifically as an `Array{T,2}` or the `Matrix{T}` type, which are equivalent.
If you are new to Julia, you can find a description of the `Array`
[here](https://docs.julialang.org/en/v1/manual/arrays/).
The first dimension contains each atomic degree of freedom, and the second dimension
contains each atom. For example, a 3D system with two atoms would have positions:
```@repl positions
using Symbolics
@variables x1, y1, z1, x2, y2, z2
r = [x1 x2;
     y1 y2;
     z1 z2]
```

!!! info "Adding external packages"

    [Symbolics](https://github.com/JuliaSymbolics/Symbolics.jl) is a package available from
    the General registry. You will have to add it to your current environment using
    `pkg> add Symbolics` to be able to reproduce this example.
    Throughout the documentation we occasionally use external packages, if you run into
    an error you will likely have to `add` the package before being able to use it.
    Refer to the [Julia manual](https://docs.julialang.org/en/v1/stdlib/Pkg/) for further
    information on installing packages.

For a 1D system it would be necessary to create a 1x1 matrix:
```@repl positions
r = fill(x1, (1,1))
```
Velocities are handled in the same way as positions and the data structures are the same.
Usually manual initialisation like this will only be necessary for small model systems,
whereas full dimensional model systems will be read from a file instead.
This is explored in the [`Atoms` documentation](@ref atoms).

!!! tip "Ring polymer simulations?"

    We can also perform simulations using ring polymers which have multiple replicas
    of each atom, these are implemented using `Array{T,3}` where the third dimension is
    used for each ring polymer bead.
    For more information, see the ring polymer methods in the dynamics methods section.

### Models

The next ingredient required to set up the simulation is the `Model`, i.e., the potentials in which the system
evolves.
These `Model`s are provided by [NQCModels.jl](@ref), which 
is a convenient infrastructure for defining different kinds of classical and quantum models
for adiabatic and nonadiabatic dynamics.
These models can range from simple analytic potentials all the way to multi-dimensional
*ab initio* potentials.
Refer to the [NQCModels.jl](@ref) page for information on the available models
and a description of how to implement new models.
 
For now we can look at a `ClassicalModel` which provides a simple
harmonic potential energy function.

```@repl started
model = Harmonic()
```
Here, the four parameters (`m`, `ω`, `r₀` and `dofs`) for this model are shown along with
their types and default values.
These values can be modified by specifying a new value in the constructor.
For example for `m`:
```@repl started
model = Harmonic(m=0.4)
```

!!! tip "Check out Parameters.jl"

    Many of the models use [Parameters.jl](https://github.com/mauro3/Parameters.jl)
    to provide convenient keyword constructors and formatted printing for the models.
    The `Harmonic` model above is defined using the `@with_kw` macro from `Parameters.jl`
    to give it a set of default parameters.
    Each of these can be modified by specifying a new value using keyword arguments in the
    constructor as demonstrated above.

Classical models implement two functions to calculate the total energy and the forces,
respectively: `potential(model, R)` and `derivative(model, R)` along with their 
in-place counterparts, `potential(model, V, R)!` and `derivative!(model, D, R)`.

Let's try these out and take a look at the results:
```@repl started
V = potential(model, hcat(25.0))
D = derivative(model, hcat(25.0))
```

!!! note "Why hcat?"

    All models accept an `R::AbstractMatrix` for the argument representing the positions
    of the particles in the system.
    These are structured such that `size(R) = (dofs, natoms)` where `dofs` is the number
    of degrees of freedom for each atom, and `natoms` is the number of atoms in the 
    simulation.

    Since this is a 1D model, we use `hcat` to quickly create a 1x1 matrix.

```@repl started
V = hcat(0.0)
D = hcat(0.0)
potential!(model, V, hcat(25.0))
derivative!(model, D, hcat(25.0))
```

!!! note "Why does potential!() treat V as a matrix?"

    The `potential!(model, V, R)` function treats V as a matrix for two reasons. The first reason 
    is that Julia treats numbers as immutable, meaning they can't be updated in place so it wouldn't 
    make sense to define `potential!()` acting on a number. The second reason is consistency; many models,
    `QuantumModels` in particular, define their potentials in terms of Hermitian matrices, therefore 
    defining a function that expects its potential to be input as a matrix is useful for interfacing 
    with other areas of the codebase.

Now, to make sure the model is what we expect, we can plot the potential and derivative 
using a custom plotting recipe. This looks pretty harmonic to me!

```@example started
using Plots

plot(-5:0.1:5, model)
```

!!! warning

    Plotting recipes currently only exist for 1D models. For more complex models you will
    have to handle the plotting manually.

### Simulation

To control all simulation parameters in one environment, we use the `Simulation` type
which will contain both the `Atoms` and `Model`s explained above, along with any
extra information required for the simulation.

```@repl started
sim = Simulation{Classical}(Atoms(:H), model)
```

Here, we have specified that each atom has a single degree of freedom and have not
provided a simulation cell.
`Classical` is a type parameter, and specifies the dynamics method that we want to use.
Check out [Dynamics methods](@ref) to learn about the other kinds of dynamics available.

!!! note

    Technically `Simulation(atoms, model)` is equivalent to
    `Simulation{Classical}(atoms, model)` since `Classical` is the default.

### Dynamics variables

The final ingredient before we can perform our simulation is the initial positions
and velocities of our particles.
For each dynamics type, the method `DynamicsVariables` is implemented and creates
the dynamics variables for us.
For classical dynamics we must provide a `Matrix` of velocities and of positions.
These should have `size = (dofs, natoms)`, matching the arguments of the `potential`
and `derivative` functions.
Usually the initial coordinates would have some physical significance, perhaps sampled
from a relevant distribution, but here we use random numbers for simplicity.
```@repl started
v = rand(3, 3);
r = rand(3, 3);
DynamicsVariables(sim, v, r)
```

!!! note

    Since `DifferentialEquations.jl` requires `AbstractArray`s for the dynamics variables,
    we use [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl)
    which allow us to conveniently store all the required information for different types of
    dynamics.

### Bringing it all together

We have now covered all the parts necessary to perform our first classical dynamics
simulation.
Let's quickly set up our simulation parameters using what we've learned.
Here we'll have two atoms in a harmonic potential, each with a single degree of freedom.

```@repl classical
using NQCDynamics # hide

atoms = Atoms([:H, :C])
sim = Simulation{Classical}(atoms, Harmonic(ω=50.0))
z = DynamicsVariables(sim, randn(size(sim)), randn(size(sim)))

nothing # hide
```

Now, we can finally run the trajectory using the `run_dynamics` function.
This takes three positional arguments: the simulation parameters `sim`,  the time span
we want to solve for, `tspan`, and the dynamics variables `z`.
For classical dynamics we also provide a timestep `dt` since we're using the
`VelocityVerlet` algorithm by default.

!!! note "Integration algorithms"

    Each method will default to an appropriate integration algorithm though it is possible
    to specify via a keyword argument to [`run_dynamics`](@ref) if an alternative
    algorithm is preferred.
    Refer to the [dynamics documentation](@ref dynamicssimulations) for more information.

The final keyword argument `output` is used to specify the quantities we want
to save during the dynamics.
A list of the available quantities can be found [here](@ref NQCDynamics.DynamicsOutputs).

!!! tip "Output format"

    `run_dynamics` returns a `Dictionary` from `Dictionaries.jl` that has entries containing
    the time and the output quantities saved at each time step.

```@example classical
tspan = (0.0, 50.0)
solution = run_dynamics(sim, (0.0, 50.0), z;
                                   dt=0.1, output=(OutputPosition, OutputVelocity))
```

Here you can see the output containing the time steps and the output quantities we specified.
These can be accessed directly as shown here:
```@repl classical
solution[:Time]
solution[:OutputPosition]
```

As with the models, we provide custom plotting recipes to quickly visualise the results
before performing further analysis by manually accessing the fields of the solution table.
To use these recipes, simply provide the solution to the `plot` function from `Plots.jl`
and give the name of the output quantity as the second argument.
This will only work if this quantity was specified in `run_dynamics`.
```@example classical
using Plots # hide
plot(solution, :OutputPosition)
plot!(solution, :OutputVelocity)
```

### Ensemble simulations

We have shown how to perform a single trajectory, but usually we are interested in performing many
and calculating observables using statistical methods.
Running more trajectories is as simple as providing the `trajectories` keyword to `run_dynamics`,
but we'll go through this in more detail in the [Ensemble simulations section.](@ref ensembles)

### What's next?

Now that we've covered the basics of classical dynamics, we're ready to explore the
world of nonadiabatic dynamics.
All the dynamics methods follow these patterns and anything you find elsewhere in the
documentation should now seem relatively familiar.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
