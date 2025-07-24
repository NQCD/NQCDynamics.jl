```@setup logging
@info "Expanding src/ensemble_simulations.md..."
start_time = time()
```
# [Ensemble simulations](@id ensembles)

Typically we'll be interested in computing observables based upon the statistics
obtained from many trajectories.
Technically it is possible to manually run many trajectories using the single trajectory
procedure introduced in the [Getting started](@ref) section.
However, by using the methods introduced on this page it is possible to run many trajectories
at once, using parallelism and computing ensemble observables automatically.

The key function for performing ensemble simulations is [`run_dynamics`](@ref).

```@docs
run_dynamics
```

This is the same function used to perform single trajectory simulations, but by replacing the single initial condition with a distribution
and changing the number of trajectories it is possible to run an ensemble of trajectories.
The distributions are defined such that they can be sampled to provide initial conditions for each trajectory.
The [Storing and sampling distributions](@ref nqcdistributions) page details the format the distributions must take.

Internally, the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/features/ensemble/#Performing-an-Ensemble-Simulation)
ensemble infrastructure is used to handle per trajectory parallelism.
The `ensemble_algorithm` keyword takes one of the [EnsembleAlgorithms](https://diffeq.sciml.ai/stable/features/ensemble/#EnsembleAlgorithms).
To use these, you must first add `using DiffEqBase` to your script.

## Example

To demonstrate usage of [`run_dynamics`](@ref), let's investigate different ways to calculate the time-dependent population
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
This procedure is detailed in the [Storing and sampling distributions](@ref nqcdistributions) section.
```@example ensemble
using Distributions: Normal # Import the Normal distribution

k = 10 # Initial momentum = 10
v = k / atoms.masses[1] # Convert momentum to velocity
r = Normal(-8) # Create Normal distribution centred at -8 for sampling initial position
nuclear_distribution = DynamicalDistribution(v, r, (1,1)) # Combine position and velocity
electronic_distribution = PureState(2) # Place electron in the 2nd electronic state (1st excited state)
product_distribution = nuclear_distribution * electronic_distribution
nothing # hide
```
In this case, we have used a deterministic momentum of ``10`` a.u. and a
gaussian position distribution with width 1 centered at ``-8`` a.u..
The electronic variables will be sampled such that the initial population is confined
to the second state by `PureState(2)`.

The final step before running the dynamics is to decide how to output the results.
The simplest option is to use the built-in tuple format familiar from [`run_dynamics`](@ref).
```@example ensemble
ensemble = run_dynamics(sim, (0.0, 3000.0), product_distribution;
    trajectories=20, output=OutputDiabaticPopulation)
nothing # hide
```
This is equivalent to performing single trajectories in a loop and manually re-sampling the initial conditions each time.
However, here we have been able to do this more concisely, using internal mechanisms for sampling from the `product_distribution`.

The output of this function is a vector containing the output from each trajectory.
Each entry is equivalent to the output from a call to [`run_dynamics`](@ref) and 
can be plotted by iterating through `ensemble`.
```@example ensemble
using Plots

p = plot(legend=false)
for traj in ensemble
    plot!(traj[:Time], [population[2] - population[1] for population in traj[:OutputDiabaticPopulation]])
end
p
```
This plot shows the population difference between the two states for each trajectory.
To approximate the exact quantum dynamics for this model, the average over all trajectories should be computed.
Instead of manually averaging the result, we can use `reduction=MeanReduction()` or `reduction=SumReduction()`
which will reduce the data accordingly before outputting the result:
```@example ensemble
ensemble = run_dynamics(sim, (0.0, 3000.0), product_distribution;
    trajectories=20, output=OutputDiabaticPopulation, reduction=MeanReduction(), saveat=0.0:10.0:3000.0)
plot(ensemble, :OutputDiabaticPopulation)
```

!!! note

    Here we have also specified the `saveat` keyword to ensure the output is saved at the
    same points for every trajectory, otherwise the averaging will not work.
    This is necessary because we are using an integrator with adaptive timestepping that will
    save at different points for each trajectory.

This workflow can be applied for any of the quantities defined in the [`DynamicsOutputs`](@ref NQCDynamics.DynamicsOutputs) submodule.
If we want a more complex output, such as a scattering probability or a time-correlation function,
we can provide a function to the output argument as described in the
[DifferentialEquations.jl documentation](https://diffeq.sciml.ai/stable/features/ensemble/#Building-a-Problem).
The advantage of this approach is that memory can be saved by reducing the data as the trajectories accumulate,
it also allows greater flexibility when modifying the output.

!!! note "Optimising performance"

    Due to how Julia's code precompilation works, running a function for the first time will take longer than subsequent calls to it. As a result, it's advantageous to run dynamics for a single time step to 
    force the Julia compiler to precompile everything needed to propagate dynamics. 
    NQCDynamics does this automatically, and you can disable this for shorter simulations by setting `precompile_dynamics=false` in the `run_dynamics` command. 

Inside the [`Ensembles`](@ref) submodule we define a few premade functions of this sort, but here
we can demonstrate how to reformulate the previous simulation using the alternative format.
```@example ensemble
function output_function(sol, i)
    output = zeros(2,div(3000, 50) + 1)
    for (i,u) in enumerate(sol.u)
        output[:,i] .= Estimators.diabatic_population(sim, u)
    end
    return output
end

ensemble = run_dynamics(sim, (0.0, 3000.0), product_distribution;
    trajectories=20, output=output_function, reduction=MeanReduction(), saveat=50.0)
```
This function provides us the same output as above, but this gives us the flexibility to calculate any observable we want.

Throughout the documentation, ensemble simulations like this one are used to demonstrate
many of the dynamics methods.
Now that you have understood the contents of this page, all of the ensemble simulations
will appear familiar.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
