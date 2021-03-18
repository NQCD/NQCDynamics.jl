# Ensemble Simulations

Typically we'll be interested in computing observables based upon the statistics
obtained from many trajectories.
We take advantage of the `EnsembleProblem` interface from `DifferentialEquations.jl` to
help us with this.

As usual we set up our system, this time we'll be doing FSSH dynamics:
```@example ensemble
using NonadiabaticMolecularDynamics # hide

atoms = NonadiabaticMolecularDynamics.Atoms([:H])
model = Models.TullyModelOne()
sim = Simulation{FSSH}(atoms, model; DoFs=1)
nothing # hide
```

When performing ensemble simulations, we provide our initial nuclear distribution in the
`DynamicalDistribution` type.
Let's create a distribution with a deterministic momentum of ``10`` a.u. and a
gaussian position distribution with width 1 centered at ``-5`` a.u.:
```@example ensemble
using Distributions: Normal

k = 10
distribution = InitialConditions.DynamicalDistribution(k / atoms.masses[1], Normal(-5), (1,1))
nothing # hide
```
Recall that the final argument here is required to specific the size of each sample.
Here we have one atom with one degree of freedom, so the size is `(1,1)`.

Now we can run the ensemble of trajectories and visualise the surface hopping populations.
```@example ensemble
using Plots

solution = Dynamics.run_ensemble(distribution, (0.0, 3000.0), sim; trajectories=50, output=(:state))

plt = plot()
for trajectory in solution
    plot!(trajectory, :state)
end
plt
```

For ensemble simulations, the result is the same as for a single trajectory, but returned in
a vector. The solution can then be indexed to extract the information from each trajectory.
