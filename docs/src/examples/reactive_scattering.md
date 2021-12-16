# [Reactive scattering from a metal surface](@id example-h2scattering)

The current version of our Julia implementation allow us to simulate the state-to-state
vibrational de-excitation probability on reactive scattering events at metal surfaces.
To run this kind of simulations, a set of initial positions and velocities with a specific
set of ro-vibrational quantum states ``\nu`` and ``j`` have to be initially generated
(see initial condition section/ sampling methods/ semiclassical EBK quantisation).
With a specific ro-vibrational quantum state is possible to compute different properties
after molecular collision and energy transfer with the metal surface like the vibrational
de-excitation probabilities.
A full description of our reactive system also includes other basic variables definition
like type of atoms, unit cell and electronic temperature which can be properly defined
before to run any simulation.

As shown earlier in the [EBK documentation](@ref ebk-sampling) we are able to generate
a semiclassically quantised distribution for a diatomic molecule on a collision course
with a metal surface.
Here we can follow that example with the [`H2AgModel`](@ref NNInterfaces.H2AgModel)
to prepare our initial distribution.

```@example h2scatter
using NonadiabaticMolecularDynamics
using NNInterfaces
using Unitful
using NonadiabaticMolecularDynamics.InitialConditions: QuantisedDiatomic

atoms = Atoms([:H, :H])
model = H2AgModel()
cell = PeriodicCell([11.1175 -5.5588 0.0; 0.0 9.628 0.0; 0.0 0.0 70.3079])
sim = Simulation(atoms, model; cell=cell)

ν, J = 2, 0
nsamples = 20

configurations = QuantisedDiatomic.generate_configurations(sim, ν, J;
    samples=nsamples, translational_energy=2.5u"eV", height=10)
v = first.(configurations)
r = last.(configurations)

distribution = DynamicalDistribution(v, r, (3,2))
```

Since we are interested in the dynamics only when the molecule is close to the surface,
we can use a callback to terminate the simulation early to save us some time.
This requires defining a function that returns `true` when we want the simulation to
terminate.
This means we can set our time span relatively long since we expect most simulations to
terminate before reaching the time limit.

```@example h2scatter
using Statistics: mean
using LinearAlgebra: norm

h2distance(p) = norm(p[:,1] .- p[:,2])

function termination_condition(u, t, integrator)::Bool
    R = get_positions(u)
    zcom = au_to_ang(mean(R[3,:]))
    if zcom > 8.1
        return true
    elseif au_to_ang(h2distance(R)) > 2.5
        return true
    else
        return false
    end
end

terminate = DynamicsUtils.TerminatingCallback(termination_condition)
tspan = (0.0, 420.0u"fs")
nothing # hide
```

## MDEF with the LDFA

Now that we've set up the initial distribution and some of our simulation parameters,
we can choose which form of friction we would like use.
First, let's use the LDFA implementation provided by the
[CubeLDFAModel](@ref models-cubeldfa).
This takes a `.cube` file containing the electron density and will provide the friction
during the dynamics.
Here we have given the new model our model from above, which will provide the forces.

```@example h2scatter
using CubeLDFAModel
model = LDFAModel(model, "../assets/friction/test.cube", atoms, friction_atoms=[1,2], cell)
```

Now we can pass all the variables defined so far to the `Simulation` and run multiple
trajectories using [`run_ensemble`](@ref).
```@example h2scatter
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_ensemble(sim, tspan, distribution;
    dt=0.1u"fs", output=:position, trajectories=20, callback=terminate)
```

Now we can plot the hydrogen bond length during the dynamics and see how some
of the trajectories lead to dissociation.
```@example h2scatter
using CairoMakie

f = Figure()
ax = Axis(f[1,1], xlabel="Time /ps", ylabel="H2 bond length /bohr")

for i=1:length(ensemble)
    lines!(au_to_ps.(ensemble[i].t), au_to_ang.(h2distance.(ensemble[i].position)))
end

f
```

## MDEF with neural network friction 

Above, we used the LDFA interpretation of MDEF to perform the simulation.
However, the [`H2AgModel`](@ref NNInterfaces.H2AgModel) actually provides it's own
friction tensor trained on *ab initio* data.
This can be used by simply using the model directly, without wrapping it with the
[`LDFAModel`](@ref CubeLDFAModel.LDFAModel).

```@example h2scatter
model = H2AgModel()
sim = Simulation{MDEF}(atoms, model, cell=cell, temperature=300u"K")
ensemble = run_ensemble(sim, tspan, distribution;
    dt=0.1u"fs", output=:position, trajectories=20, callback=terminate)

f = Figure()
ax = Axis(f[1,1], xlabel="Time /ps", ylabel="H2 bond length /bohr")

for i=1:length(ensemble)
    lines!(au_to_ps.(ensemble[i].t), au_to_ang.(h2distance.(ensemble[i].position)))
end

f
```
