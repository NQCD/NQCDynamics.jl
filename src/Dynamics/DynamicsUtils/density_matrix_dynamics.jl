using StructArrays: StructArray
using ComponentArrays: ComponentVector

using ....NonadiabaticMolecularDynamics: AbstractSimulation, Simulation, RingPolymerSimulation
using ..Dynamics: EhrenfestMethods, SurfaceHoppingMethods

const DensityMatrixMethod = Union{EhrenfestMethods.AbstractEhrenfest,
                                  SurfaceHoppingMethods.SurfaceHopping}

function Dynamics.set_quantum_derivative!(dσ, v, σ, sim::AbstractSimulation{<:DensityMatrixMethod})
    V = calculate_density_propagator!(sim, v)
    commutator!(dσ, V, σ, sim.calculator.tmp_mat_complex1)
    lmul!(-im, dσ)
end

function calculate_density_propagator!(sim::Simulation{<:DensityMatrixMethod}, v)
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.nonadiabatic_coupling[I]
    end
    return V
end

function calculate_density_propagator!(sim::RingPolymerSimulation{<:DensityMatrixMethod}, v)
    V = sim.method.density_propagator
    centroid_v = get_centroid(v)

    V .= diagm(sim.calculator.centroid_eigenvalues)
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

Dynamics.get_quantum_subsystem(u::ComponentVector{T}) where {T} =
    StructArray{Complex{T}}((u.σreal, u.σimag))
