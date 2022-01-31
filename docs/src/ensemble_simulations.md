# [Ensemble simulations](@id ensembles)

Typically we'll be interested in computing observables based upon the statistics
obtained from many trajectories.
Technically it is possible to manually run many trajectories using the single trajectory
procedure introduced in the [Getting started](@ref) section.
However, by using the methods introduced on this page it is possible to run many trajectories
at once, using parallelism and computing ensemble observables automatically.

The key function for performing ensemble simulations is [`run_ensemble`](@ref).

```@docs
run_ensemble
```

From the function signature displayed above it should be possible to identify the similarities to the [`run_trajectory'](@ref) function.
The `sim` and `tspan` positional arguments are the same, but the initial [`DynamicalVariables`](@ref) have been replaced by a distribution.
These distributions are defined such that they can be sampled to provide initial conditions for each trajectory.
The [Storing and sampling distributions](@ref) page details the format this distribution must take.

The output may take two distinct forms: the tuple structure familiar from [`run_trajectory`](@ref) `output=(:position, :velocity, :hamiltonian)`
or a function as described in the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/features/ensemble/#Performing-an-Ensemble-Simulation) documentation.
These options allow access to the common output quantities defined in the [`DynamicsOutputs`](@ref) module
along with specific customised output.
Already implemented in the code are a small library of existing functions of this type, but it is possible to use any Julia function in its place.
The existing ensemble outputs can be found [here](@ref `Ensembles`).

!!! warning "Reduction keyword"

    When using a functional output, the `reduction` keyword can be used to modify how the data is reduced between trajectories.
    If using the tuple output, only `reduction=:append` is currently valid.

Internally, the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/features/ensemble/#Performing-an-Ensemble-Simulation)
ensemble infrastructure is used to handle per trajectory parallelism.
The `ensemble_algorithm` keyword takes one of the [EnsembleAlgorithms](https://diffeq.sciml.ai/stable/features/ensemble/#EnsembleAlgorithms).
To use these, you must first add `using DiffEqBase` to your script.

## Example

To demonstrate usage of [`run_ensemble`](@ref), let's investigate different ways to calculate the time-dependent population
with [FSSH](@ref fssh-dynamics).

First, we set up our system using one of Tully's simple models ([Tully1990](@cite)).
```@example ensemble
using NQCDynamics

atoms = Atoms(2000)
model = TullyModelOne()
sim = Simulation{FSSH}(atoms, model)
nothing # hide
```

As mentioned above, before running the ensemble, we must prepare a distribution to generate
initial conditions for each trajectory.
This procedure is detailed in the [Storing and sampling distributions](@ref) section.
```@example ensemble
using Distributions: Normal # Import the Normal distribution

k = 10 # Initial momentum = 10
v = k / atoms.masses[1] # Convert momentum to velocity
r = Normal(-8) # Create Normal distribution centred at -8 for sampling initial position
nuclear_distribution = DynamicalDistribution(v, r, (1,1)) # Combine position and velocity
electronic_distribution = SingleState(2) # Create nonequilibrium electronic distribution
product_distribution = nuclear_distribution * electronic_distribution
nothing # hide
```
In this case, we have used a deterministic momentum of ``10`` a.u. and a
gaussian position distribution with width 1 centered at ``-8`` a.u..
The electronic variables will be sampled such that the initial population is confined
to the second state by `SingleState(2)`.

The final step before running the dynamics is to decide how to output the results.
The simplest option is to use the built-in tuple format familiar from [`run_trajectory`](@ref).
```@example ensemble
ensemble = run_ensemble(sim, (0.0, 3000.0), product_distribution;
    trajectories=20, output=:population)
nothing # hide
```
This is equivalent to performing single trajectories in a loop and manually re-sampling the initial conditions each time.
However, here we have been able to do this more concisely, using internal mechanisms for sampling from the `product_distribution`.

The output of this function is a vector containing the output from each trajectory.
Each entry is equivalent to the output from a call to [`run_trajectory`](@ref) and 
can be plotted by iterating through `ensemble`.
```@example ensemble
using Plots

p = plot(legend=false)
for traj in ensemble
    plot!(traj.t, [population[2] - population[1] for population in traj.population])
end
p
```
This plot shows the population difference between the two states for each trajectory.
To approximate the exact quantum dynamics for this model, the average over all trajectories should be computed.
Instead of manually averaging the result, we can use `reduction=:mean` or `reduction=:sum`
which will reduce the data accordingly before outputting the result:
```@example ensemble
ensemble = run_ensemble(sim, (0.0, 3000.0), product_distribution;
    trajectories=20, output=:population, reduction=:mean, saveat=0.0:10.0:3000.0)
plot(ensemble, :population)
```

!!! note

    Here we have also specified the `saveat` keyword.

This workflow can be applied for any of the quantities defined in the [`DynamicsOutputs`](@ref) submodule.
If we want a more complex output, such as a scattering probability or a time-correlation function,
we can provide a function to the output argument as described in the
[DifferentialEquations.jl documentation](https://diffeq.sciml.ai/stable/features/ensemble/#Building-a-Problem).
The advantage of this approach is that memory can be saved by reducing the data as the trajectories accumulate,
it also allows greater flexibility when modifying the output.

Inside the [`Ensembles`](@ref) submodule we define a few premade functions of this sort, but here
we can demonstrate how to reformulate the previous simulation using the alternative format.
```
function output_function(sol, i)
    output = zeros(2,length(sol.u))
    for (i,u) in enumerate(sol.u)
        output[:,i] = Estimators.diabatic_population(sim, u)
    end
    return (output, false)
end

ensemble = run_ensemble(sim, (0.0, 3000.0), product_distribution;
    trajectories=20, output=output_function, reduction=:mean)
```
This function provides us the same output as above, but here we have defined it
in a way compatible with the [DifferentialEquations.jl format](https://diffeq.sciml.ai/stable/features/ensemble/#Building-a-Problem).

The ensemble interface allows for specific outputs and reductions
that apply to these outputs.
Here, we choose to output the populations of each diabatic state and reduce by averaging
the results over all trajectories.
```@example ensemble
output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
nothing # hide
```

Now we can run the ensemble of trajectories and visualise the surface hopping populations.
The variable `saveat = n` allows to save the trajectory every `n` timepoints, so in this case,
we will be saving information about the trajectory every 10.0 a.u.
```@example ensemble
using Plots

solution = run_ensemble(sim, (0.0, 3000.0), distribution; trajectories=1e3,
    output=output, reduction=:mean, saveat=10.0)

plot(0:10:3000, [p[1,1] for p in solution], legend=false)
plot!(0:10:3000, [p[2,1] for p in solution])
ylabel!("Population difference")
xlabel!("Time")
```

If instead it is preferred to output many quantities for each trajectory, this is
also possible.
Here, the output is specified in the same way as for single trajectories.
```@example ensemble
ensemble = run_ensemble(sim, (0.0, 3000.0), distribution;
    output=(:population), trajectories=50)

p = plot(legend=(false))
for e in ensemble
    plot!(e.t, [pop[2] for pop in e.population])
end
p
```
Here, we see the population of the second diabatic state for every trajectory.
This is useful for checking the simulation is providing sensible results
before scaling up and outputting only the necessary data.
