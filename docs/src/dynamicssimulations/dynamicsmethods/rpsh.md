# Ring polymer surface hopping (RPSH)

Ring polymer surface hopping was one of the early attempts to extend
[RPMD](@ref Ring polymer molecular dynamics (RPMD)) to
the realm of nonadiabatic dynamics [Shushkov2012](@cite).
On the surface, the concept is reasonably simple. Since RPMD proceeds on a single
adiabatic surface, it should be possible to directly combine the
[FSSH](@ref Fewest-switches surface hopping (FSSH)) scheme with ring
polymer dynamics to approximately include some nuclear quantum effects in the surface
hopping dynamics.
However, there are some ambiguities surrounding the exact implementation when
considering how to couple the electronic equations to the nuclear equations and how the
velocity rescaling should be implemented.

Originally, two varieties were proposed: the bead approximation and the
centroid approximation.
The centroid approximation is the simpler of the two and involves directly replacing
the classical particle in the FSSH algorithm with the ring polymer centroid.
This means that the nonadiabatic couplings evaluated at the centroid and
the centroid velocity are used to propagate the electronic equations, and the
kinetic energy is conserved on the centroid level.
This is the version that is implemented here.

The bead approximation involves evaluating the nonadiabatic couplings for every bead
and using these contributions from every bead to propagate the electronics.
This version acts to conserve the kinetic energy for the entire ring polymer.
For further information on the specifics of the implementations, refer to
[Shushkov2012](@cite) and [Shakib2017](@cite).

## Example

In this example we can apply RPSH to the [`ThreeStateMorse`](@ref) model as in the
supporting info of [Shakib2017](@cite).

```@example rpsh
using NonadiabaticMolecularDynamics
using Unitful

atoms = Atoms(20000)
model = ThreeStateMorse()
sim = RingPolymerSimulation{FSSH}(atoms, model, 4; temperature=300u"K")
nothing # hide
```

For our initial conditions let's use a position distribution centred at 2.1 a.u.
with Boltzmann velocities at 300 K.
```@example rpsh
using Distributions: Normal

position = Normal(2.1, 1 / sqrt(20000 * 0.005))
velocity = InitialConditions.BoltzmannVelocityDistribution(300u"K" * nbeads(sim), masses(sim))
distribution = InitialConditions.DynamicalDistribution(velocity, position, size(sim); state=1)
nothing # hide
```

Now let's run an ensemble of trajectories that sample from this distribution.
For the output we shall receive the diabatic population at intervals of `t=50`
and it will be averaged over all trajectories by the `MeanReduction`.
```@example rpsh
solution = Ensembles.run_ensemble(sim, (0.0, 3000.0), distribution;
    saveat=50, trajectories=5e2,
    output=Ensembles.OutputDiabaticPopulation(sim), reduction=Ensembles.MeanReduction())
```

!!! note

    In the examples section at the end of the documentation we will return to this model
    and compare the performance of multiple methods.

Here we can see a plot that closely resembles the literature reference. The small
discrepancy that occurs at around `t=2000` is due to our use of a different way
to calculate the diabatic populations and better results in general would be
obtained if using the correct ring polymer distribution.
```@example rpsh
using Plots

plot(0:50:3000, [p[1] for p in solution.u], label="State 1")
plot!(0:50:3000, [p[2] for p in solution.u], label="State 2")
plot!(0:50:3000, [p[3] for p in solution.u], label="State 3")
xlabel!("Time /a.u.")
ylabel!("Population")
```

