# Getting started

To get started with the package we can identify the necessary ingredients to
perform a simple classical dynamics simulation and walk through how
to set up the simulation.

```@setup started
using NonadiabaticMolecularDynamics
```

### Atoms

First, we must define the particles in the simulation.
For this purpose we provide the `Atoms` type which will contain
the symbols, atomic numbers and masses for our atoms.
Technically these need not be actual atoms and be a generic particle.

If using real atoms, then they can be constructed using the chemical symbols
as a `Vector` of Julia's `Symbol` types, a `Vector{Symbol}`:
```@repl started
Atoms([:H, :C])
```
You can see that this contains two atoms labelled by their atomic numbers
with their masses in atomic units.

!!! note "Atomic units"

    Internally atomic units are used for all quantities. This makes things extra
    simple when performing nonadiabatic dynamics.

Alternatively, if not using real atoms, `Atoms` can be created using
a `Vector{<:Real}` where the provided numbers are the masses of the particles.
```@repl started
Atoms([1, 2, 3, 4, 5, 6])
```

### Representing atomic positions and velocities 

This package chooses to separate the dynamical variables from the static atomic parameters
included in the [`Atoms`](@ref NonadiabaticDynamicsBase.Atoms) type.
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
@variables x1, y1, z1, x2, z2, z3
r = [x1 x2;
     y1 y2;
     z1 z2]
```

For a 1D system it would be necessary to create a 1x1 matrix:
```@repl positions
r = fill(x1, (1,1))
```
Velocities are handled in the same way as positions and the data structures are the same.
Usually manual initialisation like this will only be necessary for small model systems,
whereas full dimensional model system will be read from a file instead.
This is explored in the [`Atoms` documentation](@ref Reading and writing atomic structures).

!!! tip "Ring polymer simulations?"

    We can also perform simulations using ring polymers which have multiple replicas
    of each atom, these are implemented using `Array{T,3}` where the third dimension is
    used for each ring polymer bead.
    For more information, see the ring polymer methods in the dynamics methods section.

### Models

The next ingredient required to set up the simulation is the `Model`.
These are how we tell the simulation about the potentials in which the system
evolves.
These `Model`s are provided by [NonadiabaticModels.jl](@ref), which 
is a convenient infrastructure for defining different kinds of models
for adiabatic and nonadiabatic dynamics.
But for now we can look at a simple `AdiabaticModel` which provides a simple
harmonic potential energy function.

```@repl started
model = Harmonic()
```
Here, _m_ _\omega_, r_0 and _dofs_ are the parameters and the default value of 
harmonic potential energy function. The value of the parameters can be changed,
for example for _m_:
```@repl started
model = Harmonic(m=0.4)
```


!!! tip "Check out Parameters.jl"

    Many of the models use [Parameters.jl](https://github.com/mauro3/Parameters.jl)
    to provide convenient keyword constructors and nice printing for the models.

Adiabatic models implement two functions to calculate the total energy and the forces,
respectively: `potential(model, R)` and `derivative(model, R)`.

Let's try these out and take a look at the results:
```@repl started
potential(model, hcat(25.0))
derivative(model, hcat(25.0))
```

!!! note "Why hcat?"

    All models accept an `R::AbstractMatrix` for the argument representing the positions
    of the particles in the system.
    These are structured such that `size(R) = (dofs, natoms)` where `dofs` is the number
    of degrees of freedom for each atom, and `natoms` is the number of atoms in the 
    simulation.

    Since this is a 1D model, we use `hcat` to quickly create a 1x1 matrix.

To make sure the model is what we expect, we can plot the potential and derivative 
using a custom plotting recipe. This looks pretty harmonic to me!

```@example started
using Plots

plot(-5:0.1:5, model)
```

!!! warning

    Plotting recipes currently only exist for 1D models. For more complex models you will
    have to handle the plotting manually.

### Simulation

To controll all simulation parameters in one environment, we use the `Simulation` type
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

Now we have all the parts necessary to perform our first classical dynamics simulation!

Let's quickly set up our simulation parameters using what we've learned.
Here we'll have two atoms in a harmonic potential, each with a single degree of freedom.

```@repl classical
using NonadiabaticMolecularDynamics # hide

atoms = Atoms([:H, :C])
sim = Simulation{Classical}(atoms, Harmonic(ω=50.0))
z = DynamicsVariables(sim, randn(size(sim)), randn(size(sim)))

nothing # hide
```

Now, we can finally run the trajectory using the `run_trajectory` function.
This takes three positional arguments: the dynamics variable (_z_), the time span
we want to solve for ((0.0, 50.0)), and the simulation parameters (_sim_).
For classical dynamics we also provide a timestep `dt` since we're using the
`VelocityVerlet` algorithm by default.
The final keyword argument `output` is used to specify the quantities we want
to save during the dynamics.
A list of the available quantities can be found [here](@ref DynamicsOutputs).

As with the models, we provide custom plotting recipes to quickly visualise the results
before performing further analysis.

```@example classical
using Plots # hide

solution = run_trajectory(z, (0.0, 50.0), sim;
                                   dt=0.1, output=(:position, :velocity))

plot(solution, :position)
plot!(solution, :velocity)
```

!!! note

    `run_trajectory` returns a `Table` from `TypedTables.jl` that has columns containing
    the time and the output quantities saved at each time.
    By default, it outputs the value of the dynamics variables into the field `u`.

### Ensemble simulations

So we have solved a single trajectory? That's pretty cool but wouldn't it be great
if we could do a whole bunch at once? Well, fortunately we can thanks to the
`run_ensemble` function. 

### What's next?

Now that we've covered the basics of classical dynamics, we're ready to explore the
world of nonadiabatic dynamics.
All the dynamics methods follow these patterns and anything you find elsewhere in the
documentation should now seem relatively familiar.
