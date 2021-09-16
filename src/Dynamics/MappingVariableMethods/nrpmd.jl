export NRPMD
export get_population

using LinearAlgebra: tr

"""
Nonadiabatic ring polymer molecular dynamics
"""
struct NRPMD{T} <: Method
    temp_q::Vector{T}
    temp_p::Vector{T}
    function NRPMD{T}(n_states::Integer) where {T}
        new{T}(zeros(n_states), zeros(n_states))
    end
end

select_algorithm(::RingPolymerSimulation{<:NRPMD}) = MInt()

function DynamicsVariables(sim::RingPolymerSimulation{<:NRPMD}, v, r, state::Integer; type=:diabatic)
    n_states = sim.calculator.model.n_states
    n_beads = length(sim.beads)
    if type == :diabatic
        qmap = zeros(n_states, n_beads)
        pmap = zeros(n_states, n_beads)

        θ = rand(n_states, n_beads) .* 2π
        qmap .= 1 ./ sqrt.(tan.(-θ) .^ 2 .+ 1)
        qmap[state,:] .*= sqrt(3)
        pmap = qmap .* tan.(-θ)
    else
        throw("Adiabatic initialisation not implemented, `type` kwarg should be `:diabatic`.")
    end
    return ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

get_mapping_positions(u::ComponentVector) = u.qmap
get_mapping_momenta(u::ComponentVector) = u.pmap
get_mapping_positions(u::ComponentVector, i) = @view get_mapping_positions(u)[:,i]
get_mapping_momenta(u::ComponentVector, i) = @view get_mapping_momenta(u)[:,i]

function motion!(du, u, sim::RingPolymerSimulation{<:NRPMD}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        for j in range(sim.atoms)
            for k=1:sim.DoFs
                D = sim.calculator.derivative[k,j,i]
                mul!(sim.method.temp_q, D, qmap)
                mul!(sim.method.temp_p, D, pmap)
                dv[k,j,i] = dot(qmap, sim.method.temp_q)
                dv[k,j,i] += dot(pmap, sim.method.temp_p)
                dv[k,j,i] -= tr(D)
                dv[k,j,i] /= -2sim.atoms.masses[j]
            end
        end
    end
    apply_interbead_coupling!(dv, get_positions(u), sim)
end

function set_mapping_force!(du, u, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        mul!(get_mapping_positions(du, i), V, get_mapping_momenta(u, i))
        mul!(get_mapping_momenta(du, i), V, get_mapping_positions(u, i))
        get_mapping_momenta(du, i) .*= -1
    end
end

function get_diabatic_population(::RingPolymerSimulation{<:NRPMD}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    sum(qmap.^2 + pmap.^2 .- 1; dims=2) / 2size(qmap, 2)
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::RingPolymerSimulation{<:NRPMD}, u)
    r = get_positions(u)
    v = get_velocities(u)
    Calculators.evaluate_potential!(sim.calculator, r)
    V = sim.calculator.potential

    H = get_spring_energy(sim, r) + evaluate_kinetic_energy(sim.atoms.masses, v)
    for i in range(sim.beads)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        H += 0.5 * (pmap'V[i]*pmap + qmap'V[i]*qmap - tr(V[i]))
    end
    H
end
