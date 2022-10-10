"""
    Ensembles

This module provides the main function [`run_ensemble`](@ref Ensembles.run_ensemble).
This serves to run multiple trajectories for a given simulation type, sampling
from an initial distribution.
"""
module Ensembles

using SciMLBase: SciMLBase
using Dictionaries: Dictionary
using UnitfulAtomic: austrip

using NQCBase: NQCBase
using NQCDynamics:
    AbstractSimulation,
    DynamicsMethods

export run_dynamics

include("selections.jl")

include("reductions.jl")
export SumReduction
export MeanReduction
export AppendReduction
export FileReduction

"""
    EnsembleSaver{F<:Tuple}

Store a tuple of functions with the signature `f(sol)` where `sol` is a DiffEq solution object.
`EnsembleSaver` will evaluate each of these functions and return the result in a `Dictionary`.
"""
struct EnsembleSaver{F<:Tuple}
    functions::F
end

function (output::EnsembleSaver)(sol, i)
    out = Dictionary{Symbol,Any}([:Time], [sol.t])
    return (evaluate_output_functions!(out, sol, i, output.functions...), false)
end

function evaluate_output_functions!(out::Dictionary, sol, i, f::F, fs...) where {F}
    if f isa Function
        name = nameof(f)
    else
        name = nameof(typeof(f))
    end

    insert!(out, name, f(sol, i))
    return evaluate_output_functions!(out, sol, i, fs...)
end
evaluate_output_functions!(out::Dictionary, sol, i) = out

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
function run_dynamics(
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

    prob_func = Selection(distribution, selection, trajectories)

    u0 = sample_distribution(sim, distribution)
    problem = DynamicsMethods.create_problem(u0, tspan, sim)

    output_func = EnsembleSaver(output)

    ensemble_problem = SciMLBase.EnsembleProblem(problem; prob_func, output_func, reduction)

    if trajectories == 1
        @info "Performing 1 trajectory."
    else
        @info "Performing $trajectories trajectories."
    end

    stats = @timed SciMLBase.solve(ensemble_problem, algorithm, ensemble_algorithm; trajectories, kwargs...)
    @info "Finished after $(stats.time) seconds."

    if trajectories == 1
        return stats.value.u[1]
    else
        return stats.value.u
    end
end

run_ensemble(args...; kwargs...) = throw(error("""
`run_ensemble` has been replaced, use `run_dynamics` instead.
`run_dynamics` unifies single trajectory and ensemble simulations into one function.
Refer to the `run_dynamics` docstring for more information.
"""
))

end # module
