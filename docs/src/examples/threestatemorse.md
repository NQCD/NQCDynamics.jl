# Time-dependent populations with the ThreeStateMorse model

In this example we investigate the time-dependent populations of the three state
morse model parametrised to describe photodissociation processes ([Coronado2001](@cite)).
This reference introduces three versions of this model with different parameter sets.
Our `ThreeStateMorse` model matches model C from the reference.

First let's visualise the diabats and couplings for the model.
You can see two regions where the diabats cross with non-zero coupling where we can expect
to see population transfer.
```@example threestatemorse
using NQCDynamics
using CairoMakie

x = range(2, 12, length=200)
model = ThreeStateMorse()
V = potential.(model, x)

fig = Figure()
ax = Axis(fig[1,1], xlabel="Nuclear coordinate /a.u.", ylabel="Potential energy /a.u.")

lines!(ax, x, [v[1,1] for v in V], label="State 1")
lines!(ax, x, [v[2,2] for v in V], label="State 2")
lines!(ax, x, [v[3,3] for v in V], label="State 3")

lines!(ax, x, [v[1,2] for v in V], label="Coupling 12")
lines!(ax, x, [v[2,3] for v in V], label="Coupling 23")
lines!(ax, x, [v[1,3] for v in V], label="Coupling 13")

xlims!(2, 12)
ylims!(0, 0.06)
axislegend(ax)

fig 
```

To this model we can apply any of the methods capable of starting the population on a single
diabatic state and returning the population as a function of time.
Here, let's use `FSSH` and `Ehrenfest`.
We can expect the nuclear quantum effects here to be minimal since the nuclear mass is
chosen to be 20000. 
```@example threestatemorse
m = 20000
atoms = Atoms(m)
nothing # hide
```

For our initial conditions, we use the Wigner distribution for a Harmonic oscillator
centred at 2.1 with a frequency of 5e-3 at a temperature of 300 K.
This distribution is chosen to mimic a thermal ground state distribution before
photoexcitation.
```@example threestatemorse
using Distributions: Normal
using Unitful, UnitfulAtomic

ω = 5e-3
β = 1/austrip(300u"K")
position = PositionHarmonicWigner(ω, β, m; centre=2.1)
velocity = VelocityHarmonicWigner(ω, β, m)
distribution = DynamicalDistribution(velocity, position, (1,1)) * PureState(1)
nothing # hide
```

Now let's run the two simulations using Ehrenfest and FSSH.
For both simulations we use the same initial distribution and average the results
using `reduction=:mean`.
[`TimeCorrelationFunctions.PopulationCorrelationFunction`](@ref) will correlate
the intial population with the final population at each timestep.

```@example threestatemorse
sim = Simulation{FSSH}(atoms, model)
fssh_result = run_ensemble(sim, (0.0, 3000.0), distribution;
    saveat=10, trajectories=1e3,
    output=TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic()),
    reduction=:mean, dt=1.0, u_init=[zeros(3,3) for i=1:length(0:10:3000)])
sim = Simulation{Ehrenfest}(atoms, model)
ehrenfest_result = run_ensemble(sim, (0.0, 3000.0), distribution;
    saveat=10, trajectories=1e3,
    output=TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic()),
    reduction=:mean, dt=1.0, u_init=[zeros(3,3) for i=1:length(0:10:3000)])

fig = Figure()
ax = Axis(fig[1,1], xlabel="Time /a.u.", ylabel="Population")

x = 0:10:3000
lines!(ax, x, [p[1,1] for p in fssh_result], label="FSSH State 1", color=:red)
lines!(ax, x, [p[1,2] for p in fssh_result], label="FSSH State 2", color=:green)
lines!(ax, x, [p[1,3] for p in fssh_result], label="FSSH State 3", color=:blue)

lines!(ax, x, [p[1,1] for p in ehrenfest_result], label="Ehrenfest State 1", color=:red, linestyle=:dash)
lines!(ax, x, [p[1,2] for p in ehrenfest_result], label="Ehrenfest State 2", color=:green, linestyle=:dash)
lines!(ax, x, [p[1,3] for p in ehrenfest_result], label="Ehrenfest State 3", color=:blue, linestyle=:dash)
axislegend(ax)

fig
```

To reduce the build time for the documentation the results here are underconverged but
already it is clear that both of these methods come close to the exact result shown by [Coronado2001](@cite).
After performing enough trajectories to converge the population dynamics,
we would be better able to judge the effectiveness of FSSH and Ehrenfest at reproducing the exact quantum dynamics for this model.
