# [Scattering probabilities for TullyModelTwo](@id examples-tully-model-two)

In this section we can roughly reproduce the results of figure 5 from [Tully1990](@cite).
This figure presents the scattering outcomes when a particle interacts with Tully's model 2
with an increasing magnitude of incident kinetic energy.

To reproduce this figure, first, let's set up our system parameters:
```@example tullymodeltwo
using NQCDynamics

sim = Simulation{FSSH}(Atoms(2000), TullyModelTwo())
```

Each data point in the figure is obtained from an ensemble average of trajectories.
We can use our [`Ensembles`](@ref) setup to run a set of trajectories for every single
momentum value.
Firstly, we can prepare the parts that will be the same for every ensemble:
```@example tullymodeltwo
using ComponentArrays: ComponentVector

output = Ensembles.OutputStateResolvedScattering1D(sim, :adiabatic)
```
Here, we are using the
[`OutputStateResolvedScattering1D`](@ref Ensembles.OutputStateResolvedScattering1D)
along with the [`MeanReduction`](@ref Ensembles.MeanReduction) which will give us
the average scattering outcome from the entire ensemble.
Each trajectory outputs the scattering outcome along with its final, adiabatic(?) state, and the reduction
computes the average over all trajectories.

Next, we can choose how many trajectories we want to perform for each ensemble, and
choose the range of momentum values:
```@example tullymodeltwo
ntraj = 500
momenta = 9:2:50
```

Since each ensemble requires different initial conditions, we will specific the timespan
and the distribution inside the loop.
Before the loop begins, we will create an empty list to store the results, and append
to this list after every iteration.
```@example tullymodeltwo
result = []
for k=momenta
    v = k / 2000
    r = -5
    tspan = (0, 2abs(r)/v)
    distribution = DynamicalDistribution(v, -5, size(sim)) * SingleState(1, Adiabatic())

    out = run_ensemble(sim, tspan, distribution;
        trajectories=ntraj, output=output, reduction=:mean)

    push!(result, out)
end

result
```

Now we can plot our simulation results. We format this plot to match figure 3 from
[Shakib2017](@cite). We manage to reproduce the FSSH results quite accurately. With more
trajectories, our approach will converge to the same result.

```@example tullymodeltwo
using CairoMakie

f = Figure()
ax = Axis(f[1,1], xlabel="Incident momentum / a.u.", ylabel="Scattering probability")

r1 = [r.reflection[1] for r in result]
t1 = [r.transmission[1] for r in result]
t2 = [r.transmission[2] for r in result]

scatter!(ax, momenta, r1; label="R1", color=:red)
scatter!(ax, momenta, t1; label="T1", color=:green)
scatter!(ax, momenta, t2; label="T2", color=:blue)
axislegend(ax)

f
```
