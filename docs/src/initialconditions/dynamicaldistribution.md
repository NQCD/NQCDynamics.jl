# Storing and sampling distributions

In order to perform ensembles of trajectories, it is useful to have a convenient way
to generate distributions of velocities and positions which can be sampled to
initialise trajectories.

We provide the [`DynamicalDistribution`](@ref)
type which can be used to store initial velocities and positions:
```@setup distribution
using NonadiabaticMolecularDynamics
```
```@example distribution
d = DynamicalDistribution(10, 5, (3, 2))
nothing # hide
``` 
Here we have created a delta distribution with fixed velocities and positions,
the final argument specifies the size of each sample.
The `(3, 2)` case shown here would be appropriate when using 2 atoms each with 3 degrees of freedom.
```@repl distribution
rand(d)
```

However, a delta distribution is not particularly useful, fortunately,
[`DynamicalDistribution`](@ref)
is flexible and each of the first two arguments can be `Real`, `Vector` or `Sampleable`.

!!! note

    - `Real`s are used whenever the same value is desired for every sample, as above.
    - `Vector`s can be provided when sampling a provided vector of configurations.
    - `Sampleable`s are provided by `Distributions.jl` and can be used when specifying an
        analytic distribution such as the Maxwell-Boltzmann distribution for velocities.

Each of these options can be composed in any combination, let's see the case where we have
an analytic distribution of positions and a preset collection of velocities:
```@example distribution
using Distributions

velocity = [[1.0 1.0;1.0 1.0], [2.0 2.0; 2.0 2.0], [3.0 3.0; 3.0 3.0]] 
position = Normal()
d = DynamicalDistribution(velocity, position, (2, 2))
rand(d)
``` 
This has generated normally distributed positions along with one of the three velocities
we provided.

The [`DynamicalDistribution`](@ref) also accepts
keyword arguments to provide extra metadata when performing nonequilibrium simulations
confined to a single electronic state.
Visit the docstring in API section by clicking the link in the previous sentence to learn
more about the available keywords.

This type is currently limited only to nonequilibrium simulations on a single electronic state.
It would be necessary to add further fields to store electronic variables if equilibrium
configurations were desired. This is something that can be added in the future.

To learn how to generate configurations to use with this type, read on to the next sections
about the included sampling methods.
