# Storing and sampling distributions

In order to perform ensembles of trajectories, it is useful to have a convenient way
to generate distributions of velocities and positions which can be sampled to
initialise trajectories.

Demonstrated here is the creation of the `DynamicalDistribution` which can be used to
provide initial velocities and positions:
```@setup distribution
using NonadiabaticMolecularDynamics
```
```@example distribution
d = InitialConditions.DynamicalDistribution(10, 5, (3, 2))
nothing # hide
``` 
Here we have created a distribution with fixed velocities and positions,
the final argument species the size of each sample.
The `(3, 2)` case shown here would be appropriate when using 2 atoms each with 3 DoFs.
```@repl distribution
rand(d)
```

However, a delta distribution is not particularly useful, fortunately, `DynamicalDistribution`
is flexible and each of the first two arguments can be `Real`, `Vector` or `Sampleable`.

- `Real`s are used whenever the same value is desired for every sample.
- `Vector`s can be provided when sampling a provided vector of configurations.
- `Sampleable`s are provided by `Distributions.jl` and can be used when specifying an
    analytic distribution such as the Maxwell-Boltzmann distribution for velocities.

Each of these options can be composed in any combination, let's see the case where we have
an analytic distribution of positions and a preset collection of velocities:
```@example distribution
using Distributions

velocity = [[1.0 1.0;1.0 1.0], [2.0 2.0; 2.0 2.0], [3.0 3.0; 3.0 3.0]] 
position = Normal()
d = InitialConditions.DynamicalDistribution(velocity, position, (2, 2))
rand(d)
``` 
