export CMM2

struct CMM2{T} <: Method
    temp_q::Vector{T}
    temp_p::Vector{T}
    function CMM2{T}(n_states::Integer) where {T}
        new{T}(zeros(n_states), zeros(n_states))
    end
end

const MeyerMillerMappingVariables{T} = ArrayPartition{T, Tuple{Matrix{T}, Matrix{T}, Vector{T}, Vector{T}}}

function MeyerMillerMappingVariables(v::AbstractMatrix, r::AbstractMatrix, n_states, state)

    qmap = zeros(n_states)
    pmap = zeros(n_states)

    θ = rand() * 2π
    qmap[state] = cos(θ) * sqrt(2)
    pmap[state] = sin(θ) * sqrt(2)

    ArrayPartition(v, r, pmap, qmap)
end

get_mapping_positions(z::MeyerMillerMappingVariables) = z.x[4]
get_mapping_momenta(z::MeyerMillerMappingVariables) = z.x[3]

function motion!(du, u, sim::Simulation{<:CMM2}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)

    velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u::MeyerMillerMappingVariables, sim::Simulation{<:CMM2})

    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    qmap = get_mapping_positions(u)
	pmap = get_mapping_momenta(u)
    for I in eachindex(dv)
        D = sim.calculator.derivative[I]
        mul!(sim.method.temp_q, D, qmap)
        mul!(sim.method.temp_p, D, pmap)
        dv[I] = dot(qmap, sim.method.temp_q)
        dv[I] += dot(pmap, sim.method.temp_p)
    end
    lmul!(-1/2, dv)
    divide_by_mass!(dv, sim.atoms.masses)
end

function set_mapping_force!(du, u::MeyerMillerMappingVariables, sim::Simulation{<:CMM2})
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    V = sim.calculator.potential
    mul!(get_mapping_positions(du), V, get_mapping_momenta(u))
    mul!(get_mapping_momenta(du), V, get_mapping_positions(u))
    lmul!(-1, get_mapping_momenta(du))
end

function get_diabatic_population(::Simulation{<:CMM2}, u)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    (qmap.^2 .+ pmap.^2) / 2
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::Simulation{<:CMM2}, u::MeyerMillerMappingVariables)
    r = get_positions(u)
    v = get_velocities(u)

    kinetic = evaluate_kinetic_energy(sim.atoms.masses, v)

    Calculators.evaluate_potential!(sim.calculator, r)
    V = sim.calculator.potential
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    potential = (pmap'V*pmap + qmap'V*qmap) / 2

    H = kinetic + potential
    return H
end
