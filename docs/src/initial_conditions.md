# Initial Conditions

In order to perform ensembles of trajectories, it is useful to have a convenient way
to generate distributions of velocities and positions which can be sampled to
initialise trajectories.

## Storing and sampling the distributions
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

## Generating thermal distributions
Often it is desirable to begin with a Boltzmann distributed set of velocities and positions.
Within NonadiabaticMolecularDynamics, these can be obtained using either Monte Carlo sampling,
or Langevin dynamics.

### Metropolis-Hastings Monte Carlo
This section will demonstrate how to obtain a thermal distribution in a simple
model system.

```@setup monte
using NonadiabaticMolecularDynamics
using Plots
```
First we set up the system in the usual way, here we're using an NO molecule with
a harmonic interaction between the atoms.
Notice that we use `Unitful.jl` to specify the temperature.
```@example monte
using Unitful

atoms = Atoms([:N, :O])
model = Models.DiatomicHarmonic(1.0)

sim = Simulation{Classical}(atoms, model; temperature=300u"K")
nothing # hide
```

Then we have to specify the parameters for the Monte Carlo simulation and perform the sampling.
`Δ` contains the step sizes for each of the species, `R0` the initial geometry and `passes` the
number of monte carlo passes we perform (`passes*n_atoms` steps total).
```@example monte
Δ = Dict([(:N, 0.1), (:O, 0.1)])
R0 = [1.0 0.0; 0.0 0.0; 0.0 0.0]
passes = 1000
output = MetropolisHastings.run_monte_carlo_sampling(sim, R0, Δ, passes)
nothing # hide
```

Output has three fields: the acceptance rates for each species and the energies and geometries
obtained during sampling.
```@repl monte
output.acceptance
```
```@example monte
plot(output.energy)
xlabel!("Step") # hide
ylabel!("Energy") # hide
```

We can calculate the distance between each atom and plot the bond length throughout the sampling.
```@example monte
using LinearAlgebra
plot([norm(R[:,1] .- R[:,2]) for R in output.R])
xlabel!("Step") # hide
ylabel!("Bond length") # hide
```

The result of this simulation seamlessly interfaces with the `DynamicalDistribution`
presented in the previously section and `output.R` can be readily passed to provide
the position distribution.
The Monte Carlo sampling does not include velocities but these can be readily
obtained from the Maxwell-Boltzmann distribution.
