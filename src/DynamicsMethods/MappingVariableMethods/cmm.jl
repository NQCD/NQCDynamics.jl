
struct CMM2{T} <: DynamicsMethods.Method
    temp_q::Vector{T}
    temp_p::Vector{T}
    function CMM2{T}(n_states::Integer) where {T}
        new{T}(zeros(n_states), zeros(n_states))
    end
end

function NonadiabaticMolecularDynamics.Simulation{CMM2}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.Simulation(atoms, model, CMM2{T}(NonadiabaticModels.nstates(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:CMM2}, v, r, states::Integer; type=:diabatic)

    qmap = zeros(n_states)
    pmap = zeros(n_states)

    θ = rand() * 2π
    qmap[state] = cos(θ) * sqrt(2)
    pmap[state] = sin(θ) * sqrt(2)

    ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

function motion!(du, u, sim::Simulation{<:CMM2}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u, sim::Simulation{<:CMM2})

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

function set_mapping_force!(du, u, sim::Simulation{<:CMM2})
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

function DyanmicsUtils.classical_hamiltonian(sim::Simulation{<:CMM2}, u)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)

    kinetic = DynamicsUtils.classical_kinetic_energy(sim.atoms.masses, v)

    Calculators.evaluate_potential!(sim.calculator, r)
    V = sim.calculator.potential
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    potential = (pmap'V*pmap + qmap'V*qmap) / 2

    H = kinetic + potential
    return H
end
