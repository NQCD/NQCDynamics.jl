
struct SpinMappingW{T} <: DynamicsMethods.Method
    γ::T
    R²::T
    function SpinMappingW{T}(n_states::Integer) where {T}
        γ = 2/n_states * (sqrt(n_states + 1) - 1)
        R² = 2*sqrt(n_states + 1)
        new{T}(γ, R²)
    end
end

function NQCDynamics.Simulation{SpinMappingW}(atoms::Atoms{T}, model::Model; kwargs...) where {T}
    NQCDynamics.Simulation(atoms, model, SpinMappingW{T}(NQCModels.nstates(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:SpinMappingW}, v, r, electronic::PureState{Diabatic})
    F = NQCModels.nstates(sim)

    γ = sim.method.γ
    θ = rand(F) .* 2π
    X = cos.(θ)
    P = sin.(θ)
    for i=1:F
        if i == electronic.state
            radius = sqrt(2 + γ)
        else
            radius = sqrt(γ)
        end
        X[i] *= radius
        P[i] *= radius
    end
    return ComponentVector(v=v, r=r, qmap=P, pmap=X)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:SpinMappingW}, v, r, electronic::AbstractVector{<:Integer})
    F = NQCModels.nstates(sim)

    γ = sim.method.γ
    θ = rand(F) .* 2π
    X = cos.(θ)
    P = sin.(θ)
    for i=1:F
        if i ∈ electronic
            radius = sqrt(2 + γ)
        else
            radius = sqrt(γ)
        end
        X[i] *= radius
        P[i] *= radius
    end
    return ComponentVector(v=v, r=r, qmap=P, pmap=X)
end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:SpinMappingW}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
    return nothing
end

function acceleration!(dv, u, sim::Simulation{<:SpinMappingW})

    r = DynamicsUtils.get_positions(u)

    # Set state-independent force
    NQCModels.state_independent_derivative!(sim.cache.model, dv, r)
    lmul!(-1, dv)

    NQCDynamics.NQCCalculators.update_cache!(sim.cache, r) # Ensure derivative is updated
    ∂V = NQCCalculators.get_derivative(sim.cache, r)
    X = get_mapping_positions(u)
    P = get_mapping_momenta(u)
    γ = sim.method.γ

    # Set state-dependent force
    for I in CartesianIndices(dv)
        for i in eachindex(X,P)
            dv[I] -= ∂V[I][i,i] * (X[i]^2 + P[i]^2 - γ) / 2
            for j=i+1:length(X)
                dv[I] -= ∂V[I][j,i] * (X[i]*X[j] + P[i]*P[j])
            end
        end
    end

    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses) # Convert force to acceleration
    return nothing
end

function set_mapping_force!(du, u, sim::Simulation{<:SpinMappingW})
    NQCDynamics.NQCCalculators.update_cache!(sim.cache, NQCDynamics.get_positions(u))
    r = DynamicsUtils.get_positions(u)
    V = NQCCalculators.get_potential(sim.cache, r)
    mul!(get_mapping_positions(du), V, get_mapping_momenta(u))
    mul!(get_mapping_momenta(du), V, get_mapping_positions(u))
    lmul!(-1, get_mapping_momenta(du))
    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:SpinMappingW}, u::ComponentVector)
    r = DynamicsUtils.get_positions(u)
    X = get_mapping_positions(u)
    P = get_mapping_momenta(u)
    NQCDynamics.NQCCalculators.update_cache!(sim.cache, r) # Ensure potential is updated
    V = NQCCalculators.get_potential(sim.cache, r)
    γ = sim.method.γ

    potential = NQCModels.state_independent_potential(sim.cache.model, r)
    for i in eachindex(X,P)
        potential += V[i,i] * (X[i]^2 + P[i]^2 - γ) / 2
        for j=i+1:length(X)
            potential += V[j,i] * (X[i]*X[j] + P[i]*P[j])
        end
    end
    return potential
end

function Estimators.diabatic_population(sim::Simulation{<:SpinMappingW}, u)
    X = get_mapping_positions(u)
    P = get_mapping_momenta(u)
    γ = sim.method.γ

    population = zero(X)
    for i in eachindex(X, P)
        population[i] = (X[i]^2 + P[i]^2 - γ) / 2
    end
    return population
end
