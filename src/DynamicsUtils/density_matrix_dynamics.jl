using StructArrays: StructArray
using ComponentArrays: ComponentArrays
using LinearAlgebra: diagm, mul!
using NQCDynamics: RingPolymers
using NQCModels: nstates
using RingPolymerArrays: get_centroid
using NQCDistributions: ElectronicDistribution, Adiabatic, Diabatic, density_matrix, FermiDiracState
using NQCCalculators

function set_quantum_derivative! end

function calculate_density_matrix_propagator!(propagator, v, d, eigenvalues)
    fill!(propagator, zero(eltype(propagator)))
    for (i, I) in enumerate(diagind(propagator))
        propagator[I] = eigenvalues[i]
    end

    for I in eachindex(v)
        @. propagator -= im * v[I] * d[I]
    end
    return propagator
end

"""
    commutator!(C, A, B)

Calculate C = AB - BA.
"""
function commutator!(C, A, B)
    mul!(C, B, A) # C = BA
    mul!(C, A, B, 1, -1) # C = AB - C
    return nothing
end

get_quantum_subsystem(u::ComponentArrays.ComponentVector{T}) where {T} = StructArray{Complex{T}}((u.σreal, u.σimag))

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Diabatic},
    cache::Abstract_QuantumModel_Cache,
    r
)

    diabatic_density = density_matrix(electronics, nstates(cache))
    return transform_density!(diabatic_density, cache, r, :to_adiabatic)
end

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Adiabatic},
    cache::Abstract_QuantumModel_Cache,
    r
)

    if electronics isa FermiDiracState
        eigen = NQCCalculators.get_eigen(cache, r)
        adiabatic_density = density_matrix(electronics, eigen.values)
    else
        adiabatic_density = density_matrix(electronics, nstates(cache))
    end

    return adiabatic_density
end

function transform_density!(
    density::AbstractMatrix, cache::Abstract_QuantumModel_Cache, r, direction
)
    U = evaluate_transformation(cache, r)
    if direction === :to_diabatic
        U = U'
    elseif !(direction === :to_adiabatic)
        throw(ArgumentError("`direction` $direction not recognised."))
    end
    density .= U' * density * U
    return density
end

function evaluate_transformation(cache::Abstract_QuantumModel_Cache, r)
    eigs =  NQCCalculators.get_eigen(cache, r)
    return eigs.vectors
end

function evaluate_transformation(cache::RingPolymer_QuantumModel_Cache, r)
    centroid_eigs = NQCCalculators.get_centroid_eigen(cache, r)
    return centroid_eigs.vectors
end

function evaluate_transformation(cache::RingPolymer_QuantumFrictionModel_Cache, r)
    centroid_eigs = NQCCalculators.get_centroid_eigen(cache, r)
    return centroid_eigs.vectors
end
