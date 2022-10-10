# NQCDistributions.jl

# Storing and sampling distributions

In order to perform ensembles of trajectories, it is useful to have a convenient way
to generate distributions of velocities and positions which can be sampled to
initialise trajectories.
The [NQCDistributions.jl](https://github.com/NQCD/NQCDistributions.jl) package contains the types and functions that seek
to address this requirement as painlessly as possible. 

For quantum classical nonadiabatic dynamics simulations, the initial distributions contain
both nuclear and electronic degrees of freedom.

!!! note

    Currently, we allow for *product distributions* only, where the nuclear and electronic distributions are separable.
    In the future it would be great to remove this restriction, if you are interested, please open an issue on GitHub.

This page describes the types that can be used to represent nuclear and electronic distributions
and demonstrates how they can be combined into a product distribution.

## Nuclear Distributions

### DynamicalDistribution

When handling distributions for the nuclear degrees of freedom,
the [`DynamicalDistribution`](@ref) type can be used to store initial velocities and positions:
```@setup distribution
using NQCDynamics
```
```@example distribution
d = DynamicalDistribution(10, 5, (3, 2))
nothing # hide
``` 
Here, we have created a delta distribution with fixed velocities and positions,
the final argument specifies the size of each sample.
The `(3, 2)` case shown here would be appropriate when using 2 atoms each with 3 degrees of freedom.
```@repl distribution
rand(d)
```

[`DynamicalDistribution`](@ref) is flexible and each of the first two arguments can be `Real`, `Vector` or `Sampleable`.

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

### VelocityBoltzmann

When performing equilibrium simulations it is often desirable to initialise trajectories
when thermal velocities.
These can be obtained for each atom from a gaussian distribution of the appropriate
width, or alternatively, using the [`VelocityBoltzmann`](@ref) distribution which simplifies
the process.
This takes the temperature, masses and size of the system and ensures the samples you
obtain are of the correct shape:
```@example boltzmannvelocity
using NQCDynamics
using Unitful

velocity = VelocityBoltzmann(300u"K", rand(10), (3, 10))
rand(velocity)
```
This can be handed directly to the [`DynamicalDistribution`](@ref) when Boltzmann
velocities are required.
```@example boltzmannvelocity
distribution = DynamicalDistribution(velocity, 1, (3, 10))
rand(distribution)
```

### Wigner distributions
For harmonic oscillator systems, we have implemented the analytic Wigner distribution.
These are just mormal distributions of the appropriate width but can be accessed easily
as in the following:
```@repl wigner
using NQCDistributions 

omega = 1.0;
beta = 1e-3;
mass = 10;

dist = PositionHarmonicWigner(omega, beta, mass, centre=0.0)
rand(dist)
dist = VelocityHarmonicWigner(omega, beta, mass)
rand(dist)
```
These can also be given to the [`DynamicalDistribution`](@ref) since they are just
univariate normal distributions.

## Electronic distributions

For nonadiabatic dynamics, the initial electronic variables must also be sampled.
For this, we can use an [`ElectronicDistribution`](@ref NQCDistributions.ElectronicDistribution)
which will tell our simulation how we want to sample the initial variables.
Currently, two of these are provided, the [`PureState`](@ref) and the [`MixedState`](@ref).
The [`PureState`](@ref) is used for nonequilibrium simulations when the population
is confined to a single state, whereas [`MixedState`](@ref) allows for a mixed state
distribution.

```@repl electronicdistribution
using NQCDistributions 

PureState(1, Diabatic())
PureState(2, Adiabatic())
MixedState([1, 2], Diabatic())
```

These structs contain only the minimal information about the distributions, whereas the sampling
of the distribution is handled separately by each of the different methods.

## Sampling the nuclear distribution

To learn how to generate configurations to use with the [`DynamicalDistribution`](@ref),
read on to the next sections about the included sampling methods.
