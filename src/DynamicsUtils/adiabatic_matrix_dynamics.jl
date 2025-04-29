using StructArrays: StructArray
using ComponentArrays: ComponentArrays
using LinearAlgebra: diagm, mul!
using NQCDynamics: RingPolymers
using NQCModels: nstates
using RingPolymerArrays: get_centroid
using NQCDistributions: ElectronicDistribution, Adiabatic, Diabatic, density_matrix, FermiDiracState
using ..Calculators: AbstractDiabaticCalculator, DiabaticCalculator, LargeDiabaticCalculator, RingPolymerDiabaticCalculator

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

get_quantum_subsystem(u::ComponentArrays.ComponentVector{T}) where {T} =
    StructArray{Complex{T}}((u.σreal, u.σimag))

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Diabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    diabatic_density = density_matrix(electronics, nstates(calculator))
    return transform_density!(diabatic_density, calculator, r, :to_adiabatic)
end

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Adiabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    if electronics isa FermiDiracState
        eigen = Calculators.get_eigen(calculator, r)
        adiabatic_density = density_matrix(electronics, eigen.values)
    else
        adiabatic_density = density_matrix(electronics, nstates(calculator))
    end

    return adiabatic_density
end

function transform_density!(
    density::AbstractMatrix, calculator::AbstractDiabaticCalculator, r, direction
)
    U = evaluate_transformation(calculator, r)
    if direction === :to_diabatic
        U = U'
    elseif !(direction === :to_adiabatic)
        throw(ArgumentError("`direction` $direction not recognised."))
    end
    density .= U' * density * U
    return density
end

function evaluate_transformation(calculator::Union{DiabaticCalculator, LargeDiabaticCalculator}, r)
    eigs =  Calculators.get_eigen(calculator, r)
    return eigs.vectors
end

function evaluate_transformation(calculator::RingPolymerDiabaticCalculator, r)
    centroid_eigs = Calculators.get_centroid_eigen(calculator, r)
    return centroid_eigs.vectors
end
