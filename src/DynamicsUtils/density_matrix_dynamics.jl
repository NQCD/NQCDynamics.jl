using StructArrays: StructArray
using ComponentArrays: ComponentArrays
using LinearAlgebra: diagm, mul!
using NQCDynamics: RingPolymers

function set_quantum_derivative! end

function calculate_density_matrix_propagator!(sim::Simulation, v)
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigen.values)
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.nonadiabatic_coupling[I]
    end
    return V
end

function calculate_density_matrix_propagator!(sim::RingPolymerSimulation, v)
    V = sim.method.density_propagator
    centroid_v = RingPolymers.get_centroid(v)

    V .= diagm(sim.calculator.centroid_eigen.values)
    for I in eachindex(centroid_v)
        @. V -= im * centroid_v[I] * sim.calculator.centroid_nonadiabatic_coupling[I]
    end
    return V
end

function commutator!(out, A, B, tmp)
    mul!(out, A, B)
    mul!(tmp, B, A)
    out .-= tmp
    return nothing
end

get_quantum_subsystem(u::ComponentArrays.ComponentVector{T}) where {T} =
    StructArray{Complex{T}}((u.σreal, u.σimag))
