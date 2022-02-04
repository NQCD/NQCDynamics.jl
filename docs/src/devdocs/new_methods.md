
# Contributing a new method

A key goal of NQCDynamics.jl is to provide an accessible toolkit for implementing new nonadiabatic dynamics methods.
This page details the steps you must take in order to create a new dynamics method.
The existing methods are stored inside the [`DynamicsMethods`](@ref) submodule, and this is where new methods should be implemented.
Technically, it is possible to implement new methods completely separately from the package by importing and extending
the relevant functions but if you would like to include your method in the package, it should be added within this submodule.

!!! note

    Generally, each method has its own file, though similar methods are grouped into submodules and share functionality across files.
    For example, multiple surface hopping methods have been implemented in the submodule [`DynamicsMethods.SurfaceHoppingMethods`](@ref),
    where some functions are shared across the files. 

## Basic implementation

To implement a new dynamics method, the necessary steps are:

### Create a new subtype of `DynamicsMethods.Method`.

This `Method` acts as an extra parameter inside the simulation that allows us to specify
any extra information needed for our dynamics method.
This can be a good place to include any temporary arrays and parameters for the simulation.
Refer to the [Julia manual section on Composite Types](https://docs.julialang.org/en/v1/manual/types/#Composite-Types)
to learn how this types are created.
Here, our type is called `MyMethod` and we have included an `a` parameter that will influence our dynamics:
```@example mymethod
using NQCDynamics

struct MyMethod <: DynamicsMethods.Method
    a::Float64
end
```

### Implement `DynamicsVariables`

This function returns an [`AbstractArray`](https://docs.julialang.org/en/v1/manual/arrays/)
of the variables to be used as the initial condition for the simulation.
The array should contain all of the variables that will change during the dynamics.

!!! tip "Choosing an array format"
    The only constraint on the array type is that they are [`AbstractArray`](https://docs.julialang.org/en/v1/manual/arrays/)s.
    It could be a simple matrix or vector, but usually we use [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl) to structure the variables.
    The `ComponentVector` allows us to collect variables of different types into a convenient format to perform dynamics.

For classical dynamics, this would include only the positions and velocities,
however, for [FSSH](@ref fssh-dynamics) we must also include the continuous electronic variables
and the discrete state.

!!! note "Discrete variables"
    Some methods such as [FSSH](@ref fssh-dynamics) have discontinuous variables, like the current occupied state.
    Discrete variables be handled separately using [DEDataArrays.jl](https://github.com/SciML/DEDataArrays.jl).
    For surface hopping methods, we have the `SurfaceHoppingVariables` type that uses this to combine a `ComponentVector`
    containing the continuous variables and the discrete state label.

For our new method, `MyMethod`, we implement the `DynamicsVariables` function and return a `ComponentVector` containing
the velocities, positions and extra variables `x`.
Inside this function we are free to take any inputs and manipulate them before returning the result.
As an example, suppose that our `x` variables are randomly generated each time we run the dynamics,
this could be done as follows:
```@example mymethod
using ComponentArrays

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:MyMethod}, v, r, k)
    return ComponentVector(v=v, r=r, x=rand()*k)
end
```
Here, we take the velocities `v`, positions `r` and assign them to the output as we would for classical dynamics,
but we also generate a random between 0 and `k`, where `k` was given as input.

### Implement `motion!(du, u, sim, t)`

This function should fill `du` with the time-derivative of the dynamics variables `u` in the
usual way expected by [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/tutorials/ode_example/#Example-2:-Solving-Systems-of-Equations).
We use the in-place version, where each element of `du` is filled with the time derivative of
the correponding element in `u`.

Inside the [`DynamicsUtils`](@ref) submodule there are some useful functions like [`velocity!`](@ref DynamicsUtils.velocity!) and
[`divide_by_mass!`](@ref DynamicsUtils.divide_by_mass!) which can handle some of the common parts of the `motion!` function.
You are free to perform whatever manipulations you like inside this function, but note that
`motion!` is a performance critical function, called numerous times during the simulation,
so you should attempt to minimise allocations inside this function.

!!! note
    By convention in Julia, functions that end with the `!` modify at least one of their arguments.

```@example mymethod
function DynamicsMethods.motion!(du, u, sim::Simulation{MyMethod}, t)

    DynamicsUtils.velocity!(du.r, u.v, u.r, sim, t) # Set du.r equal to the velocity

    # Set the acceleration of the particles
    du.v .= -sim.method.a .* u.r # Use the `a` parameter we stored in the `method`.
    DynamicsUtils.divide_by_mass!(du.v, sim.atoms.masses) # Divide du.v by the mass

    du.x .= 1 ./ u.x # Set time derivative of `x`.

    return nothing # The return of this function is not used so the return is unimportant
end
```

Here we have set the time derivative of the positions equal to the velocity,
the time derivative of the velocities equal to the acceleration where the force
involves the parameter `a`.
Finally, the time derivative of the extra `x` variable is also set.

### Solve a trajectory

To perform a simulation with our new method, we can write a script in the usual format
and run the dynamics.
In this script, we have a single atom with a mass of 1 with a single degree of freedom.
We match this by initialising the positions and velocities equal to random 1x1 matrices. 
The `a` parameter of the method has been set equal to 2.0, and the initial
value of `x` has been set to 0.5.

```@example mymethod
sim = Simulation(Atoms(1), Free(), MyMethod(2.0))
u = DynamicsVariables(sim, rand(1,1), rand(1,1), [0.5])

sol = run_trajectory(u, (0.0, 10.0), sim, output=(:position, :velocity, :u))
```

!!! note

    In the definition of our `motion!` method, we have accessed only the `atoms` field
    of the simulation.
    This means that the `model` we pass to the `Simulation` constructor is not used.
    Generally the `model` is accessed through the calculator interface and examples
    of its usage can be found by referring the implementations of the existing methods.

To visualise the result we can plot each of the quantities from the output table:
```@example mymethod
using Plots

plot(sol, :position, label="Position")
plot!(sol, :velocity, label="Velocity")
plot!(sol, :u, label="u", legend=true)
ylabel!("Value(t)")
```

!!! note 

    The additional `x` parameter that we created cannot be accessed in the output tuple
    by name as with `position` and `velocity` since it is not a standard quantity.
    Instead, we request `u` which contains all of the dynamical variables.
    In the plot, two of the lines labelled `u` overlap the position and velocity result.
    The unique line labelled `u` is the `x` variable.
    When implementing your method, if you want to add new output quantities you should do
    this inside the [`DynamicsOutputs`](@ref) submodule.

## Advanced tips

### Is there a custom algorithm you can implement?

Some dynamics methods have special algorithms that are tailored to the specific problem
and achieve better performance than the general algorithms include in
`DifferentialEquations.jl`.
For example, ring polymer methods typically use a symplectic scheme to solve for the
internal modes of the ring polymer, allowing much larger timesteps.
The `DifferentialEquations.jl` framework provides a simple interface for adding new
algorithms, check out the [developer documentation](https://devdocs.sciml.ai/dev/)
to learn how it works.
You can also find some examples of custom algorithms in the `DynamicsMethods.IntegrationAlgorithms` module.
