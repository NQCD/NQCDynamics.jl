
"Singleton type for labelling states as diabatic."
struct Diabatic end
"Singleton type for labelling states as adiabatic."
struct Adiabatic end

"""
    ElectronicDistribution{S}

Abstract type for distributions of electronic degrees of freedom only.
"""
abstract type ElectronicDistribution{S} end

"""
    PureState{S} <: ElectronicDistribution{S}

Electronic distribution for representing a system confined to a single state.
"""
struct PureState{S} <: ElectronicDistribution{S}
    state::Int
    statetype::S
end
PureState(state) = PureState(state, Diabatic())

function density_matrix(d::PureState, nstates)
    density = zeros(nstates, nstates)
    density[d.state, d.state] = 1
    return density
end

"""
    MixedState{T,S} <: ElectronicDistribution{S}

Electronic distribution for representing a mixed state with non-zero population in multiple states.
"""
struct MixedState{T,S} <: ElectronicDistribution{S}
    populations::Vector{T}
    statetype::S
end
MixedState(state) = MixedState(state, Diabatic())

function density_matrix(d::MixedState)
    density = zeros(length(d.populations), length(d.populations))
    density[diagind(density)] .= d.populations
    return density
end

"""
    FermiDiracState{S,T,A} <: ElectronicDistribution{S}

Electronic distribution for Fermions following Fermi-Dirac distribution.
"""
struct FermiDiracState{S,T,A} <: ElectronicDistribution{S}
    fermi_level::T
    β::T
    statetype::S
    available_states::A
end
function FermiDiracState(fermi_level, temperature; statetype=Adiabatic(), available_states=Colon())
    return FermiDiracState(austrip(fermi_level), 1/austrip(temperature), statetype, available_states)
end
fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ - μ)))

function density_matrix(d::FermiDiracState, eigenvalues)
    density = zeros(length(eigenvalues), length(eigenvalues))
    density[diagind(density)] .= fermi.(eigenvalues, d.fermi_level, d.β)
    return density
end
