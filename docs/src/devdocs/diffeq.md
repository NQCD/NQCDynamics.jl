
# DifferentialEquations.jl integration

This package is built directly on top of the established
[DifferentialEquations.jl](https://diffeq.sciml.ai/dev/index.html)
that houses a vast array of features.
Due to the tight-knit integration of our dynamics, we can immediately exploit many of
these features to save us a lot of work.
This page details how we have used DifferentialEquations.jl to make our lives easier.

## Callbacks

[Callbacks](https://diffeq.sciml.ai/dev/features/callback_functions/#callbacks) allow
us to introduce extra code during the dynamics without needing to meddle with the
integration code directly.
On the developer side, [Callbacks] is the mechanism used for the saving in the
[`run_trajectory`](@ref) function and the surface hopping procedure during FSSH.
The user can also write their own callbacks and give these to any of the dynamics functions
to manipulate the progress of the dynamics or introduce their own saving mechanism.

We also provide a few pre-made callbacks which can be given to the dynamics functions.
These are the [`TerminatingCallback`](@ref DynamicsUtils.TerminatingCallback), for terminating the simulation early,
and the [`CellBoundaryCallback`](@ref DynamicsUtils.CellBoundaryCallback)
that can be used to ensure the atoms obey the periodicity of the simulation cell.

Here, we can show how these callbacks can be used in tandem to
alter the course of the simulation. Let's look at a classical dynamics simulation without any extra callbacks:
```@example callbacks
using NQCDynamics
using Plots

atoms = Atoms(:C)
model = NQCModels.Harmonic()
cell = PeriodicCell(hcat(50))
sim = Simulation(atoms, model; cell=cell)

z = DynamicsVariables(sim, hcat(1.0), zeros(1,1))

solution = run_trajectory(z, (0.0, 300), sim; dt=0.1, output=:position)
plot(solution, :position, label="No callbacks")
```

Now we can introduce callbacks and observe the difference:
```@example callbacks
solution = run_trajectory(z, (0.0, 300), sim; callback=DynamicsUtils.CellBoundaryCallback(), dt=0.1, output=:position)
plot!(solution, :position, label="Cell boundary" )

using DiffEqBase: CallbackSet
terminate(u, t, integrator) = t > 100
callbacks = CallbackSet(DynamicsUtils.CellBoundaryCallback(), DynamicsUtils.TerminatingCallback(terminate))
solution = run_trajectory(z, (0.0, 300), sim; callback=callbacks, dt=0.1, output=:position)
plot!(solution, :position, label="Cell + termination")
```
See how the callbacks have altered the dynamics? The atom no longer leaves
the simulation cell, and the termination caused the simulation to exit early. 
