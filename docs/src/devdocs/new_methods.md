
# Contributing a new method

## Basic implementation

To implement a new dynamics method a few simple steps should be followed:

### Create a new subtype of `DynamicsMethods.Method`.

This `Method` acts as an extra parameter inside the simulation that allows us to specify
any extra information needed for our dynamics method. Here we haven't added any extra
fields but this can be a good place to include any temporary arrays and parameters
for the simulation.

```@example mymethod
using NQCDynamics

struct MyMethod <: DynamicsMethods.Method end
```

### Implement `DynamicsVariables`

This function should return an `AbstractArray` of the variables to be used as the initial
condition for the simulation.

!!! tip
    `ComponentArrays` is very useful for this. They provide the `ComponentVector` which
    makes it easy to structure and access the dynamics variables.

```@example mymethod
using ComponentArrays

function DynamicsMethods.DynamicsVariables(sim::MyMethod, v, r)
    ComponentVector(v=v, r=r)
end
```

### Implement `motion!(du, u, sim, t)`

This should fill `du` with the time-derivative of the dynamics variables `u` in the
usual way expected by `DifferentialEquations.jl`.

!!! tip
    Inside the `DynamicsUtils` submodule there are some useful functions like `velocity!` and
    `divide_by_mass!` which can handle some of the common parts of the `motion!` function.

```@example mymethod
function DynamicsMethods.motion!(du, u, sim::Simulation{MyMethod}, t)
    du.r .= u.v
    du.v .= -u.r
end
```

### Solve a trajectory

With this very simple setup we can immediately perform a simple simulation.
In this form, it is just a simple wrapper around an `OrdinaryDiffEq.jl` solver, but
by using the other parts of the framework like the `Models` and the `Calculator` interface
it becomes straightforward to implement new dynamics methods.

```@example mymethod
sim = Simulation(Atoms(1), Free(), MyMethod())
u = DynamicsVariables(sim, rand(1,1), rand(1,1))

sol = run_trajectory(u, (0.0, 10.0), sim, output=(:position, :velocity))

using Plots

plot(sol, :position)
plot!(sol, :velocity)
```

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
You can see some examples of this in practice in the `DynamicsMethods.IntegrationAlgorithms` module.
