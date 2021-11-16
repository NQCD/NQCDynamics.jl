
using LinearAlgebra: tr, mul!, dot, Diagonal
using ComponentArrays: ComponentVector
using Distributions: Normal
using NonadiabaticMolecularDynamics: RingPolymers
using NonadiabaticModels: nstates
using NonadiabaticMolecularDynamics.NonadiabaticDistributions: SingleState, Diabatic

"""
    NRPMD{T} <: DynamicsMethods.Method

Nonadiabatic ring polymer molecular dynamics
Uses Meyer-Miller-Stock-Thoss mapping variables for electronic degrees of freedom and
ring polymer formalism for nuclear degrees of freedom.

```jldoctest
RingPolymerSimulation{NRPMD}(Atoms(:H), DoubleWell(), 10)

# output

RingPolymerSimulation{NRPMD{Float64}}:
 
  Atoms{1, Float64}([:H], UInt8[0x01], [1837.4715941070515])
 
  DoubleWell{Int64, Int64, Int64, Int64}
  mass: Int64 1
  ω: Int64 1
  γ: Int64 1
  Δ: Int64 1
 
  with 10 beads.
```
"""
struct NRPMD{T} <: DynamicsMethods.Method
    γ::T
    temp_q::Vector{T}
    temp_p::Vector{T}
    function NRPMD{T}(n_states::Integer, γ) where {T}
        new{T}(γ, zeros(n_states), zeros(n_states))
    end
end

function NonadiabaticMolecularDynamics.RingPolymerSimulation{NRPMD}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; γ=0.5, kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.RingPolymerSimulation(atoms, model, NRPMD{T}(NonadiabaticModels.nstates(model), γ), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::RingPolymerSimulation{<:NRPMD}, v, r, electronic::SingleState{Diabatic})
    F = NonadiabaticModels.nstates(sim)
    n_beads = length(sim.beads)

    θ = rand(F, n_beads) .* 2π
    qmap = cos.(θ)
    pmap = sin.(θ)
    for i=1:F
        if i == electronic.state
            qmap[i, :] .*= sqrt(2 + sim.method.γ)
            pmap[i, :] .*= sqrt(2 + sim.method.γ)
        else
            qmap[i, :] .*= sqrt(sim.method.γ)
            pmap[i, :] .*= sqrt(sim.method.γ)
        end
    end

    return ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

function DynamicsMethods.motion!(du, u, sim::RingPolymerSimulation{<:NRPMD}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_derivative!(sim.calculator, DynamicsUtils.get_positions(u))
    for I in CartesianIndices(dv)
        qmap = get_mapping_positions(u, I[3])
        pmap = get_mapping_momenta(u, I[3])
        D = sim.calculator.derivative[I]
        D̄ = tr(D) / nstates(sim)
        Dtraceless = D - Diagonal(fill(D̄, nstates(sim)))
        mul!(sim.method.temp_q, Dtraceless, qmap)
        mul!(sim.method.temp_p, Dtraceless, pmap)
        dv[I] = dot(qmap, sim.method.temp_q)
        dv[I] += dot(pmap, sim.method.temp_p)
        dv[I] += 2D̄
    end
    lmul!(-1/2, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)

    DynamicsUtils.apply_interbead_coupling!(dv, DynamicsUtils.get_positions(u), sim)
end

function set_mapping_force!(du, u, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        V̄ = tr(V) / nstates(sim)
        Vtraceless = V - Diagonal(fill(V̄, nstates(sim)))
        mul!(get_mapping_positions(du, i), Vtraceless, get_mapping_momenta(u, i))
        mul!(get_mapping_momenta(du, i), Vtraceless, get_mapping_positions(u, i))
        lmul!(-1, get_mapping_momenta(du, i))
    end
end

function Estimators.diabatic_population(sim::RingPolymerSimulation{<:NRPMD}, u)
    K = NonadiabaticModels.nstates(sim)
    population = zeros(K)
    for i in 1:RingPolymers.nbeads(sim)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        for j in axes(qmap, 1)
            population[j] += (qmap[j]^2 + pmap[j]^2 - sim.method.γ) / 2
        end
    end
    return population / RingPolymers.nbeads(sim)
end

function DynamicsUtils.classical_hamiltonian(sim::RingPolymerSimulation{<:NRPMD}, u)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    spring = RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, r)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, v)

    potential = zero(eltype(u))
    Calculators.evaluate_potential!(sim.calculator, r)
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        V̄ = tr(V) / nstates(sim)
        Vtraceless = V - Diagonal(fill(V̄, nstates(sim)))
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        potential += 0.5 * (pmap'Vtraceless*pmap + qmap'Vtraceless*qmap) + V̄
    end

    return kinetic + spring + potential
end
