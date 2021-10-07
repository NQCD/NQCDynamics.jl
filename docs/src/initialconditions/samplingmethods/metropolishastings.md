# [Thermal Metropolis-Hastings Monte Carlo](@id mhmc-sampling)

Metropolis-Hastings Monte Carlo is a popular method for sampling the canonical
distribution for a molecular system.
Our implementations uses [`AdvancedMH.jl`](https://github.com/TuringLang/AdvancedMH.jl)
from the [Turing](https://turing.ml/stable/) organisation.

For a classical `Simulation`, the algorithm involves proposing new configurations in a
random walk starting from an initial configuration.
These are accepted or rejected based upon the Metropolis-Hasting criteria.
The result is a Markov chain that samples the canonical distribution.

## Example

```@example mh
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions: ThermalMonteCarlo

sim = Simulation(Atoms([:H, :H, :H, :H, :H]), Harmonic(); temperature=15)
r0 = zeros(size(sim))
chain = ThermalMonteCarlo.run_advancedmh_sampling(sim, r0, 1e4, Dict(:H=>1); move_ratio=0.5)
```

```@repl mh
Estimators.@estimate potential_energy(sim, chain)
sim.temperature / 2 * 5
```

## Legacy version

Prior to the use of [`AdvancedMH.jl`](https://github.com/TuringLang/AdvancedMH.jl),
an alternative version of the algorithm was implemented that works for both classical
and ring polymer systems.
This is currently still included in the code but should be regarded as deprecated and
will likely be removed/combined with the [`AdvancedMH.jl`](https://github.com/TuringLang/AdvancedMH.jl)
version.

Here, we use the legacy version to obtain a thermal distribution in a simple
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
model = DiatomicHarmonic(1.0)

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
output = InitialConditions.MetropolisHastings.run_monte_carlo_sampling(sim, R0, Δ, passes)
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
presented in the previous section and `output.R` can be readily passed to provide
the position distribution.
The Monte Carlo sampling does not include velocities but these can be readily
obtained from the Maxwell-Boltzmann distribution.
