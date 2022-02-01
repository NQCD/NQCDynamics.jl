# [Ring polymer surface hopping (RPSH)](@id rpsh-dynamics)

Ring polymer surface hopping was one of the early attempts to extend
[RPMD](@ref rpmd-dynamics) to
the realm of nonadiabatic dynamics [Shushkov2012](@cite).
On the surface, the concept is reasonably simple. Since RPMD proceeds on a single
adiabatic surface, it should be possible to directly combine the
[FSSH](@ref fssh-dynamics) scheme with ring
polymer dynamics to approximately include nuclear quantum effects in the surface
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
[Shushkov2012](@cite) describes both the centroid and bead approximations, [Shakib2017](@cite)
uses the centroid approximation.

## Example

In this example we can apply RPSH to the [`ThreeStateMorse`](@ref) model as shown in the
supporting info of [Shakib2017](@cite). This model has a single particle with mass of 20000 a.u.
and we use 4 beads for the ring polymer.

```@example rpsh
using NQCDynamics
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
velocity = BoltzmannVelocityDistribution(300u"K" * nbeads(sim), masses(sim), size(sim))
distribution = DynamicalDistribution(velocity, position, size(sim)) * SingleState(1)
nothing # hide
```

Now let's run an ensemble of trajectories that sample from this distribution.
For the output we will receive the diabatic population at intervals of `t=50`
and it will be averaged over all trajectories by the `:mean` keyword.
```@example rpsh
solution = run_ensemble(sim, (0.0, 3000.0), distribution;
    saveat=50, trajectories=5e2, dt=1,
    output=TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic()),
    reduction=:mean, u_init=[zeros(3,3) for i=1:length(0:50:3000)])
```

!!! note

    In the examples section at the end of the documentation we will return to this model
    and compare the performance of multiple methods.

Here we plot diabatic population of each state as a function of time.
The result closely resembles the literature reference ([Shakib2017](@cite)). The small
discrepancy that occurs at around `t=2000` is due to our use of an alternative method
to calculate the diabatic populations.
A discussion on this topic is available from [landry2013](@cite).
```@example rpsh
using Plots

plot(0:50:3000, [p[1,1] for p in solution], label="State 1")
plot!(0:50:3000, [p[1,2] for p in solution], label="State 2")
plot!(0:50:3000, [p[1,3] for p in solution], label="State 3")
xlabel!("Time /a.u.")
ylabel!("Population")
```

For our simulation we are using a Normal distribution to initialise our ring polymer configuration.
Since ring polymer surface hopping has not been rigorously derived, this choice is somewhat arbitrary
and it is possible that better results could be achieved using a free ring polymer distribution instead.
[Welsch2016](@cite) provides a theoretical description of how nonequilibrium simulations using RPMD
should be performed. This techniques here should likely be applied to RPSH too.

