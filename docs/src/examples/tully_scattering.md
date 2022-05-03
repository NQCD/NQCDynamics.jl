# [Scattering probabilities for TullyModelTwo](@id examples-tully-model-two)

In this section we aim to reproduce the results of Fig. 5 from [Tully1990](@cite).
This figure presents the scattering outcomes when a particle interacts with Tully's model 2
with an increasing magnitude of incident kinetic energy.

To reproduce this figure, first, let's set up our system parameters:
```@example tullymodeltwo
using NQCDynamics

sim = Simulation{FSSH}(Atoms(2000), TullyModelTwo())
```

!!! note Atomic units

    Recall that all of are units are atomic by default, this mass of 2000 is similar to
    that of a hydrogen atom.

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
Each trajectory outputs the scattering outcome along with its final adiabatic state, and the reduction
computes the average over all trajectories.

Next, we can choose how many trajectories we want to perform for each ensemble, and
choose the range of momentum values:
```@example tullymodeltwo
ntraj = 500
momenta = 9:2:50
```

!!! note "Range notation"

    Here we uses [Julia's range operator](https://docs.julialang.org/en/v1/base/math/#Base.::)
    to generate a set of values from 9 to 50 with a spacing of 2: 9, 11, 13, ..., 49.
    The final value of 50 is not included since a step size of 2 starting from 9 allows
    us to include only odd numbers.

Since each ensemble requires different initial conditions, we will specify the trajectory timespan
and the distribution inside the loop.
Before the loop begins, we will create an empty list to store the results, and append
to this list after every iteration.
The distribution we create produces initial conditions where each trajectory has momentum `k`
and starts at a position of `-5`. 
```@example tullymodeltwo
result = []
for k in momenta # Iterate through each momentum value
    v = k / sim.atoms.masses[1] # Starting velocity
    r = -5 # Starting position
    tspan = (0, 2abs(r)/v)
    distribution = DynamicalDistribution(v, -5, size(sim)) * PureState(1, Adiabatic())

    out = run_ensemble(sim, tspan, distribution;
        trajectories=ntraj, output=output, reduction=:mean,
        u_init=ComponentVector(reflection=zeros(2), transmission=zeros(2))
    )

    push!(result, out)
end

result
```

!!! tip "Adaptive timespan"

    Since the trajectories with larger momentum will exit the scattering region sooner,
    we scale the timespan to save computational time.
    Using `tspan = (0, 2abs(r)/v)` allows enough time such that a particle will be able
    to travel a total distance of `2r` at a constant velocity of v.
    This is sufficient to ensure the particle has left the interaction region.
    Alternatively, we could [define a callback to terminate the simulation early](@ref devdocs-callbacks).

Now we can plot our simulation results. We format this plot to match Fig. 3 from
[Shakib2017](@cite) which also reproduces Fig. 5 from [Tully1990](@cite).
We manage to reproduce the FSSH results quite accurately by visual comparison, though a larger number of trajectories
would lead to better convergence.
Since all of the examples run during the documentation build, we use a minimal number
of trajectories to optimise the build time.

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

As in [Shakib2017](@cite), R1, T1, T2 refer to reflection on state 1, transmission on
state 1 and transmission on state 2 respectively.
For this model, surface hopping is successful in closely approximating the exact quantum
result, especially at higher momentum values.
Refer to [Tully1990](@cite) and [Shakib2017](@cite) for a detailed discussion of the results.
