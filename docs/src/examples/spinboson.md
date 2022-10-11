# Ohmic spin-boson nonequilibrium population dynamics

The spin-boson model is widely used as a model for condensed phase quantum dynamics.
It is defined by a system-bath Hamiltonian where the system is a 2-state spin
coupled to a bath of harmonic oscillators.
This example shows how to perform nonequilibrium population dynamics with the spin-boson
model using a bath characterised by the Ohmic spectral density.
We will be using model B from the work of [gao2020](@cite).

Our boson bath will have 100 oscillators, each with a mass of 1.
Here, we also set up the model with the ohmic density and the parameters that
match up with our reference ([gao2020](@cite)).
The ohmic density is given a cutoff frequency of 2.5 and a Kondo parameter of 0.09.
The model is symmetric, with the energy bias between states equal to 0.0, and
the coupling between states set to 1.

```@example spinboson
using NQCDynamics
N = 100
atoms = Atoms(fill(1, N))
β = 5
T = 1 / β
density = OhmicSpectralDensity(2.5, 0.09)
model = SpinBoson(density, N, 0.0, 1.0)
nothing # hide
```

## Initial conditions

For the initial conditions, we will sample directly from a Wigner distribution for
the nuclear degrees of freedom.
Since our nuclear degrees of freedom are harmonic, the Wigner distribution has an
analytic form and we can use the distributions included in the package.
The `position` and `velocity` variables we create here are `Matrix`s of `Normal` distributions,
which are shaped to match the system size `(1, N)`.
Inside the `DynamicalDistribution` they will provide samples that match the size of the system.
The initial electronic state is confined to 1 with `PureState(1)`.

```@example spinboson
position = reshape([PositionHarmonicWigner(ω, β, 1) for ω in model.ωⱼ], 1, :)
velocity = reshape([VelocityHarmonicWigner(ω, β, 1) for ω in model.ωⱼ], 1, :)
distribution = DynamicalDistribution(velocity, position, (1, 100)) * PureState(1)
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
output = TimeCorrelationFunctions.PopulationCorrelationFunction(fssh, Diabatic())
ensemble_fssh = run_dynamics(fssh, (0.0, 20.0), distribution;
    saveat=saveat, trajectories=100, output, reduction=MeanReduction())
output = TimeCorrelationFunctions.PopulationCorrelationFunction(ehrenfest, Diabatic())
ensemble_ehrenfest = run_dynamics(ehrenfest, (0.0, 20.0), distribution;
    saveat=saveat, trajectories=100, output, reduction=MeanReduction())
nothing # hide
```

Here, we can see the population difference between the two states.
```@example spinboson
using Plots
plot(saveat, [p[1,1] - p[1,2] for p in ensemble_fssh[:PopulationCorrelationFunction]], label="FSSH")
plot!(saveat, [p[1,1] - p[1,2] for p in ensemble_ehrenfest[:PopulationCorrelationFunction]], label="Ehrenfest")
xlabel!("Time /a.u.")
ylabel!("Population difference")
```

The exact result for this model, along with various mapping methods can be found 
in the work of [gao2020](@cite).
We can see that even with just 100 trajectories, our Ehrenfest result closely matches theirs.
The FSSH is quite clearly underconverged with only 100 trajectories due to the discontinuous
nature of the individual trajectories.
Feel free to try this for yourself and see what the converged FSSH result looks like!
