module Postprocess
using NQCDynamics: AbstractSimulation, Ensembles

export apply_output_functions

struct FakeProblem
    p::AbstractSimulation
end

struct FakeSolution
    t::AbstractArray
    u::AbstractArray
    p::FakeProblem
end

"""
    apply_output_functions(sol::FakeSolution, output_functions; savetime::Bool=true)

Evaluates output functions on a DifferentialEquations.jl solution object or fake solution object
generated using a defined Simulation and a `DynamicsVariables` type.

Basically equivalent to running `run_dynamics()` with the same output functions, but without
doing the dynamics simulation again.
"""
function apply_output_functions(sol, output_functions; savetime::Bool=true)
    if !isa(output_functions, Tuple)
        output_functions = (output_functions,)
    end
    output = Ensembles.EnsembleSaver(output_functions, savetime)
    return output(sol, 1)
end

"""
    apply_output_functions(sim<:AbstractSimulation, u_type::AbstractVector, t_type::AbstractVector, output_functions; savetime::Bool=true)

Evaluates output functions on a defined Simulation, a `DynamicsVariables` type output and a time-type output.

Basically equivalent to running `run_dynamics()` with the same output functions, but without
doing the dynamics simulation again.
"""
function apply_output_functions(sim :: AbstractSimulation, u_type::AbstractVector, t_type::AbstractVector, output_functions; savetime::Bool=true)
    sol = FakeSolution(t_type, u_type, FakeProblem(sim))
    return apply_output_functions(sol, output_functions; savetime=savetime)
end

end
