
using StochasticDiffEq: StochasticDiffEq
using RecursiveArrayTools: ArrayPartition
import NQCModels.FrictionModels

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
struct Langevin{gt<:AbstractFloat, st<:AbstractMatrix{<:AbstractFloat}} <: AbstractMDEF
    γ::gt
    σ::st
end

function Langevin(γ::AbstractFloat, temperature::Number, masses::AbstractVector, DoFs::Integer)
    σ = sqrt.(temperature ./ repeat(Array(masses'), DoFs, 1))
    Langevin(γ, σ)
end

# StochasticDiffEq.BAOAB is currently broken, so let's do our one instead and call it MDEF.
function NQCDynamics.Simulation{Langevin}(atoms::Atoms, model::Model; γ=1, temperature=0u"K", kwargs...)
    gamma_matrix = zeros(eltype(atoms.masses), length(atoms.masses), length(atoms.masses))
    for i in 1:length(atoms.masses)
        gamma_matrix[i, i] = γ
    end
    gamma_ondevice = similar(atoms.masses, size(gamma_matrix)...)
    gamma_ondevice .= gamma_matrix

    friction_matrix = NQCModels.FrictionModels.ConstantFriction(ndofs(model), gamma_ondevice)
    combined_model = CompositeFrictionModel(model, friction_matrix)

    NQCDynamics.Simulation(atoms, combined, MDEF(atoms.masses, ndofs(model));
               temperature=temperature, kwargs...)
end

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::Simulation{<:Langevin})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end
# DynamicsMethods.select_algorithm(sim::AbstractSimulation{<:Langevin}) = StochasticDiffEq.BAOAB(;gamma=sim.method.γ)

# function friction!(du, r, sim::AbstractSimulation{<:Langevin}, t)
#     du .= sim.method.σ
# end

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
struct ThermalLangevin{T<:Real} <: DynamicsMethods.Method
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
