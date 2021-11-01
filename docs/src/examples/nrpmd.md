# Nonadiabatic ring polymer molecular dynamics -- NRPMD

NRPMD is a method for nonadiabatic dynamics that uses the ring polymer formalism
for representing the nuclear degrees of freedom and the mapping variable formalism
for the electronic variables.

```@example nrpmd
using NonadiabaticMolecularDynamics
using Plots
nothing # hide
```

Let's apply NRPMD to the spin boson model.

Our boson bath will have 4 oscillators, each with a mass of 1.
```@example nrpmd
N = 4
atoms = Atoms(fill(1, N))
T = 1 / 16
n_beads = 5
density = OhmicSpectralDensity(0.25, 0.09)
model = SpinBoson(density, N, 0.0, 1.0)
sim = RingPolymerSimulation{NRPMD}(atoms, model, n_beads; temperature=T)
nothing # hide
```

## Initial conditions

To setup the initial conditions we can sample the thermal distribution associated with the
boson bath in isolation from the spin.
This is done most easily by defining a second model that includes only the boson
degrees of freedom.
```@example nrpmd
bath = BosonBath(density, N)
nothing # hide
```

We can perform the sampling using thermal Langevin dynamics.
Since the type of dynamics is determined by the simulation type,
we must create a modified `RingPolymerSimulation`.
```@example nrpmd
langevin_sim = RingPolymerSimulation{ThermalLangevin}(atoms, bath, n_beads;
    temperature=T)

v = RingPolymerArray(zeros(1, N, n_beads))
r = RingPolymerArray(zeros(1, N, n_beads))
z = DynamicsVariables(langevin_sim, v, r)

sol = run_trajectory(z, (0.0, 10.0), langevin_sim;
    dt=0.01, output=(:position, :velocity))
```

Let's plot the velocity and position of each of the beads to check we have reasonable
sampling.
```@example nrpmd
v_plot = plot(sol, :velocity)
r_plot = plot(sol, :position)
plot(v_plot, r_plot, layout=(2,1))
```

From our results we can create a distribution that can be sampled for the NRPMD dynamics.
```@example nrpmd
distribution = DynamicalDistribution(sol.velocity, sol.position, size(sim)) * SingleState(1)
nothing # hide
```

## NRPMD dynamics

Now that we have a distribution from which we can sample our initial conditions,
we can run the ensemble of trajectories and calculate the diabatic populations.
```@example nrpmd
output = Ensembles.OutputDiabaticPopulation(sim)
reduction = Ensembles.MeanReduction()
ensemble = Ensembles.run_ensemble(sim, (0.0, 30.0), distribution;
    saveat=0.1, trajectories=10, output=output, reduction=reduction, dt=1e-3)
nothing # hide
```

```@example nrpmd
plot(0:0.1:30, [p[1] - p[2] for p in ensemble.u])
xlabel!("Time /a.u.")
ylabel!("Population difference")
```
