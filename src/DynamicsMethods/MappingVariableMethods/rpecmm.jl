using LinearAlgebra: lmul!, norm
using Parameters: Parameters
using NQCDistributions: ElectronicDistribution
using NQCDynamics: nbeads
using NQCModels: nstates

function NQCDynamics.RingPolymerSimulation{eCMM}(atoms::Atoms{T}, model::Model, n_beads::Integer; γ=0, kwargs...) where {T}
    NQCDynamics.RingPolymerSimulation(atoms, model, eCMM{T}(nstates(model), γ), n_beads; kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::RingPolymerSimulation{<:eCMM}, v, r, ::ElectronicDistribution)

    F = nstates(sim)
    n_beads = length(sim.beads)

    radius = sqrt(2 + 2F*sim.method.γ*n_beads)
    points = generate_random_points_on_nsphere(2F*n_beads, radius)
    half = length(points) ÷ 2
    qmap = points[begin:half]
    pmap = points[half+1:end]
    qmap = reshape(qmap, (F, n_beads))
    pmap = reshape(pmap, (F, n_beads))

    return ComponentVector(v=v, r=r, pmap=pmap, qmap=qmap)
end

function acceleration!(dv, u, sim::RingPolymerSimulation{<:eCMM})

    Parameters.@unpack γ, temp_q, temp_p = sim.method

    NQCCalculators.evaluate_derivative!(sim.cache, DynamicsUtils.get_positions(u))
    for I in CartesianIndices(dv)
        qmap = get_mapping_positions(u, I[3])
        pmap = get_mapping_momenta(u, I[3])
        D = sim.cache.derivative[I]
        D̄ = tr(D) / nstates(sim)
        Dtraceless = D - Diagonal(fill(D̄, nstates(sim)))
        mul!(temp_q, Dtraceless, qmap)
        mul!(temp_p, Dtraceless, pmap)
        dv[I] = dot(qmap, temp_q)
        dv[I] += dot(pmap, temp_p)
        dv[I] += 2D̄
    end
    lmul!(-1/2, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    DynamicsUtils.apply_interbead_coupling!(dv, DynamicsUtils.get_positions(u), sim)
end

function set_mapping_force!(du, u, sim::RingPolymerSimulation{<:eCMM})
    NQCCalculators.evaluate_potential!(sim.cache, DynamicsUtils.get_positions(u))
    for i in range(sim.beads)
        V = sim.cache.potential[i]
        V̄ = tr(V) / nstates(sim)
        Vtraceless = V - Diagonal(fill(V̄, nstates(sim)))
        mul!(get_mapping_positions(du, i), Vtraceless, get_mapping_momenta(u, i))
        mul!(get_mapping_momenta(du, i), Vtraceless, get_mapping_positions(u, i))
        lmul!(-1, get_mapping_momenta(du, i))
    end
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:eCMM}, u)
    r = DynamicsUtils.get_positions(u)

    NQCCalculators.evaluate_potential!(sim.cache, r)
    potential = zero(eltype(u))
    for i=1:nbeads(sim)
        V = sim.cache.potential[i]
        V̄ = tr(V) / nstates(sim)
        Vtraceless = V - Diagonal(fill(V̄, nstates(sim)))
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        potential += 0.5 * (pmap'Vtraceless*pmap + qmap'Vtraceless*qmap) + V̄
    end

    return potential
end

function inverse_mapping_kernel(qmap, pmap, γ, N)
    F = length(qmap)
    prefactor = (1 + N*F) / (1 + N*F*γ)^2
    subtractor = (1 - γ) / (1 + N*F*γ)
    return prefactor .* (qmap.^2 .+ pmap.^2)./2 .- subtractor
end

function Estimators.diabatic_population(sim::RingPolymerSimulation{<:eCMM}, u)
    population = zeros(nstates(sim))
    for i=1:nbeads(sim)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        population .+= inverse_mapping_kernel(qmap, pmap, sim.method.γ, nbeads(sim))
    end
    return population
end

function Estimators.initial_diabatic_population(sim::RingPolymerSimulation{<:eCMM}, u)
    population = zeros(nstates(sim))
    for i=1:nbeads(sim)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        population .+= mapping_kernel(qmap, pmap, sim.method.γ)
    end
    return population
end
