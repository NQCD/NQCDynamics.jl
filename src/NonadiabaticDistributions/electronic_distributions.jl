
abstract type StateType end
struct Diabatic <: StateType end
struct Adiabatic <: StateType end

abstract type ElectronicDistribution{S<:StateType} <: NonadiabaticDistribution end

struct SingleState{S} <: ElectronicDistribution{S}
    state::Int
    statetype::S
end
SingleState(state) = SingleState(state, Diabatic())

struct ElectronicPopulation{T,S} <: ElectronicDistribution{S}
    populations::Vector{T}
    statetype::S
end
ElectronicPopulation(state) = ElectronicPopulation(state, Diabatic())
