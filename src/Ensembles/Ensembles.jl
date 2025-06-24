"""
    Ensembles

This module provides the main function [`run_dynamics`](@ref run_dynamics).
This serves to run multiple trajectories for a given simulation type, sampling
from an initial distribution.
"""
module Ensembles

using SciMLBase: SciMLBase
using Dictionaries: Dictionary
using UnitfulAtomic: austrip

using NQCBase: NQCBase
using NQCCalculators
using NQCDynamics: AbstractSimulation, DynamicsMethods

export run_dynamics

include("selections.jl")

include("reductions.jl")
export SumReduction
export MeanReduction
export AppendReduction
export FileReduction

include("run_dynamics.jl")
export run_dynamics
export log_simulation_duration

"""
    EnsembleSaver{F<:Tuple}

Store a tuple of functions with the signature `f(sol)` where `sol` is a DiffEq solution object.
`EnsembleSaver` will evaluate each of these functions and return the result in a `Dictionary`.
"""
struct EnsembleSaver{F<:Tuple}
    functions::F
    savetime::Bool
end

function (output::EnsembleSaver)(sol, i)
    if output.savetime
        out = Dictionary{Symbol,Any}([:Time], [sol.t])
    else
        out = Dictionary{Symbol,Any}()
    end
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


end # module
