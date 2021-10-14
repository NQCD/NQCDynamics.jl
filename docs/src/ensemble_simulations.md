# [Ensemble simulations](@id ensembles)

Typically we'll be interested in computing observables based upon the statistics
obtained from many trajectories.

As usual we set up our system, this time we'll be doing FSSH dynamics:
```@example ensemble
using NonadiabaticMolecularDynamics # hide

atoms = NonadiabaticMolecularDynamics.Atoms([:H])
model = TullyModelOne()
sim = Simulation{FSSH}(atoms, model)
nothing # hide
```

When performing ensemble simulations, we provide our initial nuclear distribution in the
`DynamicalDistribution` type.
Let's create a distribution with a deterministic momentum of ``10`` a.u. and a
gaussian position distribution with width 1 centered at ``-8`` a.u.:
```@example ensemble
using Distributions: Normal

k = 10
v = k / atoms.masses[1]
r = Normal(-8)
distribution = InitialConditions.DynamicalDistribution(v, r, (1,1); state=2)
nothing # hide
```
!!! note

    The keyword argument `state` used here is for nonadiabatic simulations and specifies the
    initial electronic state.

The ensemble interface allows for specific outputs and reductions
that apply to these outputs.
Here we choose to output the populations of each diabatic state and reduce by averaging
the results over all trajectories.
```@example ensemble
output = Ensembles.OutputDiabaticPopulation(sim)
reduction = Ensembles.MeanReduction()
nothing # hide
```

Now we can run the ensemble of trajectories and visualise the surface hopping populations.
```@example ensemble
using Plots

solution = Ensembles.run_ensemble(sim, (0.0, 3000.0), distribution; trajectories=1e3,
    output=output, reduction=reduction, saveat=10.0)

plot(0:10:3000, [p[1] for p in solution.u], legend=false)
plot!(0:10:3000, [p[2] for p in solution.u])
ylabel!("Population difference")
xlabel!("Time")
```

If instead it is preferred to output many quantities for each trajectory, this is
also possible.
Here the output is specified in the same way as for single trajectories.
```@example ensemble
ensemble = Ensembles.run_trajectories(sim, (0.0, 3000.0), distribution;
    output=(:population), trajectories=50)

p = plot(legend=(false))
for e in ensemble
    plot!(e.t, [pop[2] for pop in e.population])
end
p
```
Here we see the population of the second diabatic state for every trajectory.
This is useful for checking the simulation is providing sensible results
before scaling up and outputting only the necessary data.
