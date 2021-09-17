
using StochasticDiffEq: StochasticDiffEq
using RecursiveArrayTools: ArrayPartition
using Unitful: @u_str
using UnitfulAtomic: austrip

struct Langevin{T<:AbstractFloat} <: DynamicsMethods.Method
    γ::T
    σ::Matrix{T}
end

function Langevin{T}(γ, temperature, masses, DoFs) where {T}
    σ = sqrt.(temperature ./ repeat(Array(masses'), Int(DoFs), 1))
    Langevin(T(γ), T.(σ))
end

function NonadiabaticMolecularDynamics.Simulation{Langevin}(atoms::Atoms{S,T}, model::Model; γ=1, temperature=0u"K", DoFs=3, kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.Simulation(atoms, model, Langevin{T}(γ, austrip(temperature), atoms.masses, DoFs);
               temperature=temperature, DoFs=DoFs, kwargs...)
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::Simulation{<:Langevin})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end
DynamicsMethods.select_algorithm(sim::AbstractSimulation{<:Langevin}) = StochasticDiffEq.BAOAB(sim.method.γ)

function friction!(du, r, sim::AbstractSimulation{<:Langevin}, t)
    du .= sim.method.σ
end

struct ThermalLangevin{T<:Real} <: DynamicsMethods.Method
    γ::T
end

function NonadiabaticMolecularDynamics.RingPolymerSimulation{ThermalLangevin}(atoms::Atoms, model::Model, n_beads::Integer; γ=1, kwargs...)
    NonadiabaticMolecularDynamics.RingPolymerSimulation(atoms, model, ThermalLangevin(γ), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(::AbstractSimulation{<:Union{ThermalLangevin, Langevin}}, v, r)
    ArrayPartition(v, r)
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:ThermalLangevin})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end
