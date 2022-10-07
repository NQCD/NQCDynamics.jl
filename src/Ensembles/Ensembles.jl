"""
    Ensembles

This module provides the main function [`run_ensemble`](@ref Ensembles.run_ensemble).
This serves to run multiple trajectories for a given simulation type, sampling
from an initial distribution.
"""
module Ensembles

using SciMLBase: SciMLBase, EnsembleProblem, solve
using RecursiveArrayTools: ArrayPartition

using NQCBase: NQCBase
using NQCDynamics:
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsUtils,
    DynamicsMethods,
    DynamicsOutputs
using RingPolymerArrays: RingPolymerArray

export run_ensemble

include("selections.jl")
include("outputs.jl")
include("reductions.jl")

export SumReduction
export MeanReduction
export AppendReduction
export FileReduction

"""
    run_ensemble(sim::AbstractSimulation, tspan, distribution;
        selection::Union{Nothing,AbstractVector}=nothing,
        output=(sol,i)->(sol,false),
        reduction=:append
        reduction::Symbol=:append,
        ensemble_algorithm=SciMLBase.EnsembleThreads(),
        algorithm=DynamicsMethods.select_algorithm(sim),
        saveat=[],
        trajectories=1,
        kwargs...
        )

Run multiple trajectories for timespan `tspan` sampling from `distribution`.
The DifferentialEquations ensemble interface is used which allows us to specify
functions to modify the output and how it is reduced across trajectories.

# Keywords

* `selection` should be an `AbstractVector` containing the indices to sample from the `distribution`. By default, `nothing` leads to random sampling.
* `output` can be a function that transforms the DiffEq solution to an output, or a tuple of output quantities as for `run_trajectory`.
* `reduction` defines how the data is reduced across trajectories. Options are `:append`, `:mean` or `:sum`.
* `ensemble_algorithm` is the algorithm from DifferentialEquations which determines which form of parallelism is used.
* `algorithm` is the algorithm used to integrate the equations of motion.
* `saveat` mirrors the DiffEq keyword and labels the time points to save the output.
* `trajectories` is the number of trajectories to perform.

This function wraps the EnsembleProblem from DifferentialEquations and passes the `kwargs`
to the `solve` function.
"""
function run_ensemble(
    sim::AbstractSimulation, tspan, distribution;
    output,
    selection::Union{Nothing,AbstractVector}=nothing,
    reduction=AppendReduction(),
    ensemble_algorithm=SciMLBase.EnsembleSerial(),
    algorithm=DynamicsMethods.select_algorithm(sim),
    trajectories=1,
    kwargs...
)

    if !(output isa Tuple)
        output = (output,)
    end

    kwargs = NQCBase.austrip_kwargs(;kwargs...)
    trajectories = convert(Int, trajectories)
    tspan = austrip.(tspan)

    prob_func = Selection(distribution, selection)

    u0 = sample_distribution(sim, distribution)
    problem = DynamicsMethods.create_problem(u0, tspan, sim)

    output_func = DynamicsOutputs.EnsembleSaver(output)

    ensemble_problem = EnsembleProblem(problem; prob_func, output_func, reduction)

    @info "Performing $trajectories trajectories."
    stats = @timed solve(ensemble_problem, algorithm, ensemble_algorithm; trajectories, kwargs...)
    @info "Finished after $(stats.time) seconds."

    return stats.value.u
end

end # module
