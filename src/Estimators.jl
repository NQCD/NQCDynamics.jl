
"""
    Estimators

Functions for computing thermal expectation values as ensemble averages.
"""
module Estimators

using NQCDynamics:
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsUtils,
    RingPolymers,
    ndofs,
    natoms,
    nbeads,
    masses,
    get_ring_polymer_temperature
using NQCCalculators
using StatsBase: mean
using ComponentArrays: ComponentVector
using RingPolymerArrays: get_centroid

"""
    @estimate f(simulation, vector)

Evaluate `f(simulation, vector[i])` for all `i` and return the average.

Can be used for any function defined in `Estimators.jl`.
"""
macro estimate(expr)
    func = expr.args[1]
    sim = expr.args[2]
    configurations = expr.args[3]
    first_config, rest, result, config = gensym(), gensym(), gensym(), gensym()

    return esc(quote
        local ($first_config, $rest) = Iterators.peel($configurations)
        local $result = Estimators.$func($sim, $first_config)
        for $config in $rest
            $result += Estimators.$func($sim, $config)
        end
        $result / length($configurations)
    end)
end

"""
    total_energy(sim, u)
"""
function total_energy(sim::AbstractSimulation, u)
    return kinetic_energy(sim, u) + potential_energy(sim, u)
end

"""
    potential_energy(sim, u)
"""
function potential_energy(sim::AbstractSimulation, u)
    return potential_energy(sim, DynamicsUtils.get_positions(u))
end

function potential_energy(sim::Simulation, r::AbstractMatrix)
    NQCCalculators.update_potential!(sim.cache, r)
    return sim.cache.potential[1]
end

function potential_energy(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    NQCCalculators.update_potential!(sim.cache, r)
    return mean(sim.cache.potential)[1]
end

"""
    kinetic_energy(sim, u)
"""
function kinetic_energy(sim::Simulation, u)
    kinetic_energy(sim, DynamicsUtils.get_velocities(u))
end

function kinetic_energy(sim::Simulation, v::AbstractMatrix)
    DynamicsUtils.classical_kinetic_energy(sim, v)
end

function kinetic_energy(sim::RingPolymerSimulation, u)
    kinetic_energy(sim, DynamicsUtils.get_positions(u))
end

function kinetic_energy(sim::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    centroid = get_centroid(r)

    NQCCalculators.update_derivative!(sim.cache, r)

    kinetic = ndofs(sim) * natoms(sim) * get_ring_polymer_temperature(sim)

    for I in CartesianIndices(r)
        kinetic += (r[I] - centroid[I[1], I[2]]) * sim.cache.derivative[I] 
    end

    return kinetic / 2nbeads(sim)
end
    
"""
    radius_of_gyration(sim, r)
"""
function radius_of_gyration(::RingPolymerSimulation, r::AbstractArray{T,3}) where {T}
    centroid = get_centroid(r)
    deviation = (r .- centroid) .^ 2
    mean_deviation = dropdims(mean(deviation; dims=3); dims=3)
    return sqrt.(mean_deviation)
end

"""
    diabatic_population
"""
function diabatic_population end

"""
    adiabatic_population
"""
function adiabatic_population end

initial_diabatic_population(sim, u) = diabatic_population(sim, u)
initial_adiabatic_population(sim, u) = adiabatic_population(sim, u)

end # module
