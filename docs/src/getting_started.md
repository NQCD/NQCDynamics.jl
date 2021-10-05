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

### Models

The next ingredient required to set up the simulation is the `Model`.
These are how we tell the simulation about the potentials in which the particles
evolve.
These `Model`s are provided by [NonadiabaticModels.jl](@ref) which 
provides a convenient infrastructure for defining different kinds of models
for nonadiabatic dynamics.
But for now we can look at a simple `AdiabaticModel` which provides a simple
harmonic potential energy function.

```@repl started
model = Harmonic()
```

!!! tip "Check out Parameters.jl"

    Many of the models use [Parameters.jl](https://github.com/mauro3/Parameters.jl)
    to provide convenient keyword constructors and nice printing for the models.

Adiabatic models implement two functions: `potential(model, R)` and 
`derivative(model, R)`.

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

To group all simulation parameters together, we use the `Simulation` type
which will contain both the `Atoms` and `Model`s mentioned previously, along with any
extra information required for the simulation.

```@repl started
sim = Simulation{Classical}(Atoms(:H), model)
```

Here we have specified that each atom has a single degree of freedom and have not
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
This takes three positional arguments: the dynamics variable, the time span
we want to solve for, and the simulation parameters.
For classical dynamics we also provide a timestep `dt` since we're using the
`VelocityVerlet` algorithm.
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
