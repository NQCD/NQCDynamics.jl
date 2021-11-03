# Ohmic spin-boson nonequilibrium population dynamics

The spin-boson model is widely used as a model for condensed phase quantum dynamics.
It is defined by a system-bath Hamiltonian where the system is a 2-state spin
coupled to a bath of harmonic oscillators.
This example shows how to perform nonequilibrium population dynamics with this
model using a bath characterised by the Ohmic spectral density.

```@example spinboson
using NonadiabaticMolecularDynamics
using Plots
nothing # hide
```

Our boson bath will have 100 oscillators, each with a mass of 1.
Here we also set up the model with the ohmic density and an appropriate set of parameters.
```@example spinboson
N = 100
atoms = Atoms(fill(1, N))
T = 1 / 5
density = OhmicSpectralDensity(2.5, 0.09)
model = SpinBoson(density, N, 0.0, 1.0)
nothing # hide
```

## Initial conditions

For the initial conditions, we will sample directly from the Wigner distribution for
the nuclear degrees of freedom.
Since our nuclear degrees of freedom are harmonic, the Wigner distribution has an
analytic form and we can use the distributions included in the package.

```@example spinboson
position = PositionHarmonicWigner.(model.ωⱼ, β, 1)
velocity = VelocityHarmonicWigner.(model.ωⱼ, β, 1)
distribution = DynamicalDistribution(velocity, position, (1, 100))) * SingleState(1)
nothing # hide
```

## Dynamics

Now that we have a distribution from which we can sample our initial conditions,
we can run ensembles of trajectories and calculate the population correlation functions.
Let's compare the results obtained using FSSH and Ehrenfest.
```@example spinboson
fssh = Simulation{FSSH}(atoms, model)
ehrenfest = Simulation{Ehrenfest}(atoms, model)

saveat = 0:0.1:20
output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
reduction = Ensembles.MeanReduction([zeros(2, 2) for _ in saveat])
ensemble = Ensembles.run_ensemble(sim, (0.0, 20.0), distribution;
    saveat=saveat, trajectories=100, output=output, reduction=reduction)
nothing # hide
```

Here we can see the population difference between the two states.
For this set of parameters, Ehrenfest outperforms FSSH and comes close to the exact quantum
result.
```@example spinboson
plot(saveat, [p[1] - p[2] for p in ensemble.u])
xlabel!("Time /a.u.")
ylabel!("Population difference")
```
