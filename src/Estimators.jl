
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
    DynamicsUtils,
    RingPolymers,
    ndofs,
    natoms,
    masses,
    get_temperature

using StatsBase: mean
using ComponentArrays: ComponentVector

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

function potential_energy(sim::Simulation, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    sim.calculator.potential
end

function potential_energy(sim::RingPolymerSimulation, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    mean(sim.calculator.potential)
end

function kinetic_energy(sim::Simulation, u)
    v = DynamicsUtils.get_velocities(u)
    sum(masses(sim)' .* v .^2) / 2
end

function kinetic_energy(sim::RingPolymerSimulation, u)
    r = DynamicsUtils.get_positions(u)
    centroid = RingPolymers.get_centroid(r)

    Calculators.evaluate_derivative!(sim.calculator, r)

    kinetic = ndofs(sim) * natoms(sim) * get_temperature(sim)

    for I in eachindex(r)
        kinetic += (r[I] - centroid[I]) * sim.calculator.derivative[I] 
    end

    return kinetic / 2nbeads(sim)
end

function diabatic_population end
function adiabatic_population end

end # module
