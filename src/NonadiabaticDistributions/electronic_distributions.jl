using LinearAlgebra: diagind

using NQCModels: NQCModels
using NQCDynamics: Calculators, RingPolymers
using .Calculators:
    AbstractDiabaticCalculator,
    DiabaticCalculator,
    RingPolymerDiabaticCalculator

abstract type StateType end
"Singleton type for labelling states as diabatic."
struct Diabatic <: StateType end
"Singleton type for labelling states as adiabatic."
struct Adiabatic <: StateType end

"""
    ElectronicDistribution{S<:StateType} <: NonadiabaticDistribution

Abstract type for distributions of electronic degrees of freedom only.
"""
abstract type ElectronicDistribution{S<:StateType} <: NonadiabaticDistribution end

"""
    SingleState{S} <: ElectronicDistribution{S}

Electronic distribution for representing a system confined to a single state.
"""
struct SingleState{S} <: ElectronicDistribution{S}
    state::Int
    statetype::S
end
SingleState(state) = SingleState(state, Diabatic())

"""
    ElectronicPopulation{T,S} <: ElectronicDistribution{S}

Electronic distribution for representing a mixed state with non-zero population in multiple states.
"""
struct ElectronicPopulation{T,S} <: ElectronicDistribution{S}
    populations::Vector{T}
    statetype::S
end
ElectronicPopulation(state) = ElectronicPopulation(state, Diabatic())

"""
    FermiPopulation{T,S} <: ElectronicDistribution{S}

Electronic distribution for Fermions following Fermi-Dirac distribution.
"""
struct FermiPopulation{T,S} <: ElectronicDistribution{S}
    fermi_level::T
    temperature::T
    statetype::S
end
FermiPopulation(fermi_level, temperature) = ElectronicPopulation(fermi_level, temperature, Adiabatic())

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Diabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    diabatic_density = initialise_density_matrix(electronics, calculator)
    return transform_density!(diabatic_density, calculator, r, :to_adiabatic)
end

function initialise_adiabatic_density_matrix(
    electronics::ElectronicDistribution{Adiabatic},
    calculator::AbstractDiabaticCalculator,
    r
)

    adiabatic_density = initialise_density_matrix(electronics, calculator)
    return adiabatic_density
end

function initialise_density_matrix(
    electronics::ElectronicDistribution, calculator::AbstractDiabaticCalculator
)

    n = NQCModels.nstates(calculator.model)
    density = Matrix{eltype(calculator)}(undef, n, n)
    fill_density!(density, electronics)
    return density
end

function fill_density!(density::AbstractMatrix, electronics::SingleState)
    fill!(density, zero(eltype(density)))
    density[electronics.state, electronics.state] = 1
    return density
end

function fill_density!(density::AbstractMatrix, electronics::ElectronicPopulation)
    fill!(density, zero(eltype(density)))
    density[diagind(density)] .= electronics.populations
    return density
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

function evaluate_transformation(calculator::DiabaticCalculator, r)
    eigs =  Calculators.get_eigen(calculator, r)
    return eigs.vectors
end

function evaluate_transformation(calculator::RingPolymerDiabaticCalculator, r)
    centroid_eigs = Calculators.get_centroid_eigen(calculator, r)
    return centroid_eigs.vectors
end
