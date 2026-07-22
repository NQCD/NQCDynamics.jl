
using StochasticDiffEq: StochasticDiffEq
using RecursiveArrayTools: ArrayPartition
using Unitful: @u_str
using UnitfulAtomic: austrip

"""
Type for performing Langevin molecular dynamics.

```jldoctest
using Unitful
sim = Simulation{Langevin}(Atoms(:H), Free(); γ=2.5, temperature=100u"K")

# output

Simulation{Langevin{Float64}}:
  Atoms{Float64}([:H], [1], [1837.4715941070515])
  Free(1)
```
"""

abstract type LangevinMethods <: DynamicsMethods.Method end

struct Langevin{T<:AbstractFloat} <: LangevinMethods
    γ::T
    σ::Matrix{T}
end

function Langevin{T}(γ, temperature, masses, DoFs) where {T}
    σ = sqrt.(temperature ./ repeat(Array(masses'), Int(DoFs), 1))
    Langevin(T(γ), T.(σ))
end

function NQCDynamics.Simulation{Langevin}(atoms::Atoms{T}, model::Model; γ=1, temperature=0u"K", kwargs...) where {T}
    NQCDynamics.Simulation(atoms, model, Langevin{T}(γ, austrip(temperature), atoms.masses, ndofs(model));
               temperature=temperature, kwargs...)
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::Simulation{<:Langevin})
    StochasticDiffEq.SDEProblem(acceleration!, friction!, u0, tspan, sim; noise_rate_prototype = zero(DynamicsUtils.get_velocities(u0)))
end

function friction!(du, r, sim::AbstractSimulation{<:Langevin}, t)
    du .= sim.method.σ
end

"""
Type for performing Langevin ring polymer molecular dynamics.

Currently there are separate types for classical and ring polymer versions of Langevin
dynamics but they should be combined. The reason they are not at the moment is that they use
different integration algorithms and require slightly different fields.

```jldoctest
using Unitful
RingPolymerSimulation{ThermalLangevin}(Atoms(:H), Free(), 10; γ=0.1, temperature=25u"K")

# output

RingPolymerSimulation{ThermalLangevin{Float64}}:
 
  Atoms{Float64}([:H], [1], [1837.4715941070515])
 
  Free(1)
  with 10 beads.
```
"""
struct ThermalLangevin{T<:Real} <: LangevinMethods
    γ::T
end

function NQCDynamics.RingPolymerSimulation{ThermalLangevin}(atoms::Atoms, model::Model, n_beads::Integer; γ=1, kwargs...)
    NQCDynamics.RingPolymerSimulation(atoms, model, ThermalLangevin(γ), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(::AbstractSimulation{<:Union{ThermalLangevin, Langevin}}, v, r)
    ArrayPartition(v, r)
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:ThermalLangevin})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end
