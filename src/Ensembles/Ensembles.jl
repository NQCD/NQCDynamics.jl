"""
    Ensembles

This module provides the main function [`run_ensemble`](@ref Ensembles.run_ensemble).
This serves to run multiple trajectories for a given simulation type, sampling
from an initial distribution.
"""
module Ensembles

using SciMLBase: SciMLBase, EnsembleProblem, solve
using RecursiveArrayTools: ArrayPartition
using TypedTables: TypedTables

using NQCBase: NQCBase
using NQCDynamics:
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsUtils,
    DynamicsMethods,
    RingPolymers,
    NonadiabaticDistributions,
    DynamicsOutputs
using NQCDynamics.NonadiabaticDistributions:
    NuclearDistribution,
    CombinedDistribution

export run_ensemble

function sample_distribution(sim::AbstractSimulation, distribution::NuclearDistribution, i)
    u = getindex(distribution, i)
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r)
end

function sample_distribution(sim::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}, distribution::NuclearDistribution, i)
    u = getindex(distribution, i)
    DynamicsMethods.DynamicsVariables(sim, RingPolymers.RingPolymerArray(u.v), RingPolymers.RingPolymerArray(u.r))
end

function sample_distribution(sim::AbstractSimulation, distribution::CombinedDistribution, i)
    u = getindex(distribution.nuclear, i)
    DynamicsMethods.DynamicsVariables(sim, u.v, u.r, distribution.electronic)
end

include("selections.jl")
include("outputs.jl")
include("reductions.jl")

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
    selection::Union{Nothing,AbstractVector}=nothing,
    output=SciMLBase.DEFAULT_OUTPUT_FUNC,
    reduction::Symbol=:append,
    ensemble_algorithm=SciMLBase.EnsembleThreads(),
    algorithm=DynamicsMethods.select_algorithm(sim),
    saveat=[],
    trajectories=1,
    u_init=nothing,
    kwargs...
)

    kwargs = NQCBase.austrip_kwargs(;kwargs...)
    trajectories = convert(Int, trajectories)
    saveat = austrip.(saveat)
    tspan = austrip.(tspan)
    (output isa Symbol) && (output = (output,))
    reduction = Reduction(reduction)

    ensemble_problem = EnsembleProblem(sim, tspan, distribution, selection, output, reduction, u_init, saveat, trajectories, kwargs)

    @info "Performing $trajectories trajectories."
    stats = @timed solve(ensemble_problem, algorithm, ensemble_algorithm; saveat, trajectories, kwargs...)
    @info "Finished after $(stats.time) seconds."

    return extract_output(stats, output, reduction, trajectories)
end

function SciMLBase.EnsembleProblem(sim::AbstractSimulation, tspan, distribution, selection, output, reduction, u_init, saveat, trajectories, kwargs)

    prob_func = Selection(distribution, selection)

    problem = create_problem(sim, distribution, tspan)
    output_func = Output(output, sim)
    if output_func isa DynamicsOutputs.EnsembleSaver
        u_init = get_u_init(reduction, saveat, kwargs, tspan, problem.u0, output_func)
    end

    return EnsembleProblem(problem; prob_func, output_func, reduction, u_init)
end

Output(output::Tuple, sim) = DynamicsOutputs.EnsembleSaver(output, sim)
Output(output, _) = output

function create_problem(sim::AbstractSimulation, distribution, tspan)
    u0 = sample_distribution(sim, distribution, 1)
    return DynamicsMethods.create_problem(u0, tspan, sim)
end

function extract_output(stats, output, reduction, trajectories)
    sol = stats.value

    if output isa Tuple
        process_output(reduction, sol, output, trajectories)
    else
        return sol.u
    end
end

function process_output(::AppendReduction, sol, output, trajectories)
    return [generate_table(output, s) for s in sol.u]
end

function process_output(::SumReduction, sol, output, trajectories)
    for i in eachindex(sol.u)
        sol.u[i][1] /= trajectories
    end
    return generate_table(output, sol.u)
end

function process_output(::MeanReduction, sol, output, trajectories)
    return generate_table(output, sol.u)
end

function generate_table(output, u)
    output = (:t, output...)
    out = [(;zip(output, quantity)...) for quantity in u]
    return TypedTables.Table(out)
end

end # module
