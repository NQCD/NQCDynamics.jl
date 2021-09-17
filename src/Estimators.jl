
"""
    Estimators

Functions for computing thermal expectation values as ensemble averages.
"""
module Estimators

using NonadiabaticMolecularDynamics:
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    Calculators,
    DynamicsUtils

using StatsBase: mean
using ComponentArrays: ComponentVector

# function potential_estimator(sim::AbstractSimulation, u::Vector)
#     V = zero(eltype(u))
#     for frame in u
#         V += evaluate_potential_energy(sim, get_positions(frame))
#     end
#     return V / length(u)
# end
# function estimate(sim::AbstractSimulation, configurations::Vector{<:ComponentVector{T}}, func::Symbol) where {T}
#     estimation = zero(eltype(configurations[1]))
#     for config in configurations
#         estimation += eval(func)(sim, config)::T
#     end
#     return estimation
# end

"""
    @estimate f(simulation, vector)

Evaluate `f(simulation, vector[i])` for all `i` and return the average.

Can be used for any function defined in `Estimators.jl`.
"""
macro estimate(expr)
    func = expr.args[1]
    sim = expr.args[2]
    configurations = expr.args[3]
    result, config = gensym(), gensym()

    return esc(quote
        local $result = 0
        for $config in $configurations
            $result += Estimators.$func($sim, $config)
        end
        $result / length($configurations)
    end)
end

function potential(sim::Simulation, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    sim.calculator.potential
end

function potential(sim::RingPolymerSimulation, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    mean(sim.calculator.potential)
end

function diabatic_population end
function adiabatic_population end

end # module
