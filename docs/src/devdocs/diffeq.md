```@setup logging
@info "Expanding src/devdocs/diffeq.md..."
start_time = time()
```

# DifferentialEquations.jl integration

NQCDynamics.jl is built directly on top of the established
[DifferentialEquations.jl](https://diffeq.sciml.ai/dev/index.html)
that provides a vast array of features.
By using DifferentialEquations.jl to perform the dynamics,
we can immediately exploit many of these features to save us a lot of work.
This page details some of the features from DifferentialEquations.jl that we have used.

## [Callbacks](@id devdocs-callbacks)

[Callbacks](https://diffeq.sciml.ai/dev/features/callback_functions/#callbacks) allow
us to introduce extra code during the dynamics without needing to meddle with the
integration code directly.
On the developer side, [Callbacks] is the mechanism used for the saving in the
[`run_dynamics`](@ref) function and the surface hopping procedure during FSSH.
The user can also write their own callbacks and give these to any of the dynamics functions
to manipulate the progress of the dynamics or introduce their own saving mechanism.

We also provide two pre-made callbacks which can be given to the dynamics functions.
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

solution = run_dynamics(sim, (0.0, 300), z; dt=1.0, output=OutputPosition)
plot(solution, :OutputPosition, label="No callbacks", legend=true)
```

Now we can introduce callbacks and observe the difference:
```@example callbacks
solution = run_dynamics(sim, (0.0, 300), z; callback=DynamicsUtils.CellBoundaryCallback(), dt=1.0, output=OutputPosition)
plot!(solution, :OutputPosition, label="Cell boundary" )

using DiffEqBase: CallbackSet
terminate(u, t, integrator) = t > 100
callbacks = CallbackSet(DynamicsUtils.CellBoundaryCallback(), DynamicsUtils.TerminatingCallback(terminate))
solution = run_dynamics(sim, (0.0, 300), z; callback=callbacks, dt=1.0, output=OutputPosition)
plot!(solution, :OutputPosition, label="Cell + termination")
```
See how the callbacks have altered the dynamics? The atom no longer leaves
the simulation cell, and the termination caused the simulation to exit early. 

The callback setup we're using is exactly that provided by DifferentialEquations.jl,
if you want more details on callbacks, please refer to their [documentation](https://diffeq.sciml.ai/dev/features/callback_functions/#callbacks).
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
