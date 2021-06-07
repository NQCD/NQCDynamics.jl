# Fewest switches surface hopping -- FSSH

FSSH models excited state dynamics by propagating the nuclei
on a single adiabatic surface at a time, switching to different surfaces
stochastically.
The stochastic switching is based upon the electronic populations obtained
by integrating the electronic Schrodinger equation alongside the nuclear dynamics.

As usual, first we import the packages.
```@example fssh
using NonadiabaticMolecularDynamics
using Plots
using Unitful
using Distributions
nothing # hide
```

To illustrate the fewest-switches dynamics we can use a popular three state Morse potential.
```@example fssh
atoms = Atoms(20000)
model = ThreeStateMorse()
sim = Simulation{FSSH}(atoms, model; DoFs=1)
nothing # hide
```

For our initial conditions let's use a position distribution centred at 2.1 a.u.
with zero velocity.
```@example fssh
position = Normal(2.1, 1 / sqrt(20000 * 0.005))
distribution = InitialConditions.DynamicalDistribution(0.0, position, (1,1); state=1)
nothing # hide
```

Now let's run an ensemble of trajectories that sample from this distribution.
```@example fssh
solution = Ensembles.run_ensemble(sim, (0.0, 3000.0), Ensembles.RandomSelection(distribution);
    saveat=50, trajectories=5e2,
    output=Ensembles.OutputDiabaticPopulation(sim), reduction=Ensembles.MeanReduction())

plot(0:50:3000, [p[1] for p in solution.u], label="State 1")
plot!(0:50:3000, [p[2] for p in solution.u], label="State 2")
plot!(0:50:3000, [p[3] for p in solution.u], label="State 3")
xlabel!("Time /a.u.")
ylabel!("Population")
```
