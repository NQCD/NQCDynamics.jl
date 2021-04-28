# Molecular dynamics with electronic friction -- MDEF

MDEF is a method designed to include the effect of low energy electron-hole pair
excitations that occur when a molecule encounters a metal.
It can be seen as a direct extension to Born-Oppenheimer molecular dynamics
and manifests itself in the form of Langevin dynamics.

First let's import the packages we need.
```@example mdef
using NonadiabaticMolecularDynamics
using Plots
using Unitful
nothing # hide
```

Here we shall model a single hydrogen atom in a harmonic potential,
where the electronic temperature is 300 K.
The friction provided by this model is a random positive definite tensor.
```@example mdef
atoms = Atoms([:H])
model = Models.FrictionHarmonic()
sim = Simulation{MDEF}(atoms, model; temperature=300u"K")
nothing # hide
```

Next we should set up our initial conditions
```@example mdef
velocity = zeros(3, 1)
position = zeros(3, 1)
z = ClassicalDynamicals(velocity, position)
nothing # hide
```

Finally we can run a single trajectory and visualise the total energy as a function of time.
```@example mdef
solution = Dynamics.run_trajectory(z, (0.0, 1u"ps"), sim, dt=0.1u"fs",
    output=(:energy, :position, :velocity))
plot(solution, :energy)
```

Now let's see what happens if we make the electronic temperature a function of time.
```@example mdef
temperature_function(t) = exp(-(t - 0.5u"ps")^2 * 20u"ps^-2") * 3000u"K"
nothing # hide
```
!!! warning

    The time argument enters this function as a `Unitful.jl` quantity, so it
    is important to make sure the unit of the return value is temperature.

Now we can re-simulate, replacing the fixed temperature with the function we have defined.
```@example mdef
sim = Simulation{MDEF}(atoms, model; temperature=temperature_function)
solution = Dynamics.run_trajectory(z, (0.0, 1u"ps"), sim, dt=0.1u"fs",
    output=(:energy, :position, :velocity))
plot(solution, :energy)
```
This time we see a peak in the energy in the middle of the simulation which coincides
with the peak in temperature at 0.5 ps.

