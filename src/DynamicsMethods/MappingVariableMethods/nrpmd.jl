
using LinearAlgebra: tr, mul!, dot
using ComponentArrays: ComponentVector
using NonadiabaticMolecularDynamics: RingPolymers

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
    temp_q::Vector{T}
    temp_p::Vector{T}
    function NRPMD{T}(n_states::Integer) where {T}
        new{T}(zeros(n_states), zeros(n_states))
    end
end

function NonadiabaticMolecularDynamics.RingPolymerSimulation{NRPMD}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.RingPolymerSimulation(atoms, model, NRPMD{T}(NonadiabaticModels.nstates(model)), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::RingPolymerSimulation{<:NRPMD}, v, r, electronic::NonadiabaticDistributions.SingleState)
    n_states = NonadiabaticModels.nstates(sim.calculator.model)
    n_beads = length(sim.beads)
    state = electronic.state
    if electronic.statetype === NonadiabaticDistributions.Diabatic()
        qmap = zeros(n_states, n_beads)
        pmap = zeros(n_states, n_beads)

        θ = rand(n_states, n_beads) .* 2π
        qmap .= 1 ./ sqrt.(tan.(-θ) .^ 2 .+ 1)
        qmap[state,:] .*= sqrt(3)
        pmap = qmap .* tan.(-θ)
    else
        throw("Adiabatic initialisation not implemented.")
    end
    return ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

get_mapping_positions(u::ComponentVector) = u.qmap
get_mapping_momenta(u::ComponentVector) = u.pmap
get_mapping_positions(u::ComponentVector, i) = @view get_mapping_positions(u)[:,i]
get_mapping_momenta(u::ComponentVector, i) = @view get_mapping_momenta(u)[:,i]

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
        mul!(sim.method.temp_q, D, qmap)
        mul!(sim.method.temp_p, D, pmap)
        dv[I] = dot(qmap, sim.method.temp_q)
        dv[I] += dot(pmap, sim.method.temp_p)
        dv[I] -= tr(D)
    end
    lmul!(-1/2, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)

    DynamicsUtils.apply_interbead_coupling!(dv, DynamicsUtils.get_positions(u), sim)
end

function set_mapping_force!(du, u, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        mul!(get_mapping_positions(du, i), V, get_mapping_momenta(u, i))
        mul!(get_mapping_momenta(du, i), V, get_mapping_positions(u, i))
        lmul!(-1, get_mapping_momenta(du, i))
    end
end

function Estimators.diabatic_population(::RingPolymerSimulation{<:NRPMD}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    sum(qmap.^2 + pmap.^2 .- 1; dims=2) / 2size(qmap, 2)
end

function DynamicsUtils.classical_hamiltonian(sim::RingPolymerSimulation{<:NRPMD}, u)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    spring = RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, r)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, v)

    potential = zero(eltype(u))
    Calculators.evaluate_potential!(sim.calculator, r)
    V = sim.calculator.potential
    for i in range(sim.beads)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        potential += 0.5 * (pmap'V[i]*pmap + qmap'V[i]*qmap - tr(V[i]))
    end

    return kinetic + spring + potential
end
