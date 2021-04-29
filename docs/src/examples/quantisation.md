# Semiclassical quantisation of a diatomic molecule

In surface science it is often of interest to investigate the change in the quantum
state of a diatomic molecule after scattering off a surface.
We can investigate this with molecular dynamics by using semiclassical quantisation
procedures to extract the quantum numbers from classical dynamics.

Firstly, we can generate a distribution initialised in a specific quantum state
```@example quantise
using NonadiabaticMolecularDynamics
using Plots

v = 2 # Vibrational quantum number
J = 0 # Rotational quantum number

atoms = Atoms([:C, :O])
model = Models.DiatomicHarmonic()
sim = Simulation(atoms, model)

configs = QuantisedDiatomic.generate_configurations(sim, v, J;
    samples=1000, translational_energy=10)

v = first.(configs)
r = last.(configs)
distribution = InitialConditions.DynamicalDistribution(v, r, (3,2))
nothing # hide
```

```@example quantise
selection = Ensembles.RandomSelection(distribution)
output = Ensembles.OutputQuantisedDiatomic(sim)
ensemble = Ensembles.run_ensemble(sim, (0.0, 100.0), selection;
    output=output, trajectories=10, dt=0.1)
(first.(ensemble.u), last.(ensemble.u))
```
Here we see that the quantum numbers do not change from their initial values through the
dynamics. This is to be expected since we are doing classical dynamics with a harmonic
potential.
This example should be improved using a better model.
