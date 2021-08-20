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

    (n_states, n_beads) = size(qmap)
    population = zeros(n_states)
    for i=1:n_beads
        for j=1:n_states
            population[j] += qmap[j,i]^2 + pmap[j,i]^2 - 1
        end
    end
    population ./= 2n_beads
    return population
end

function get_adiabatic_population(sim::RingPolymerSimulation{<:NRPMD}, u)

    Calculators.update_electronics!(sim.calculator, get_positions(u))

    U = sim.calculator.eigenvectors

    diabatic = zeros(sim.calculator.model.n_states)
    adiabatic = zero(diabatic)
    for i=1:length(sim.beads)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        Usquare = U[i].^2
        diabatic .= (qmap .^2 .+ pmap .^2 .- 1) ./ 2
        adiabatic .+= Usquare' * diabatic
    end
    adiabatic ./= length(sim.beads)

    return adiabatic
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
