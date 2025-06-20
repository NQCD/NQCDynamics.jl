using NQCBase
using NQCDynamics: nbeads, TemperatureSetting
using RingPolymerArrays: eachbead

function NQCDynamics.RingPolymerSimulation{DiabaticMDEF}(atoms::Atoms{T}, model::Model, n_beads::Integer; friction_method::FrictionEvaluationMethod, 
    temperature=0u"K", cell::AbstractCell=InfiniteCell(), solver::Symbol=:exact, kwargs...
    ) where {T}

    cache = Create_Cache(model, length(atoms), n_beads, T; friction_method=friction_method)

    # If a thermostat is provided, check it covers the whole system.
    isa(temperature, TemperatureSetting{Vector{Int}}) ? throw(DomainError(temperature, "TemperatureSetting must apply to all atoms.")) : nothing

    # If multiple TemperatureSettings are provided, check that each atom only has one thermostat applied to it.
    if isa(temperature, Vector{<:TemperatureSetting})
        indices = vcat([thermostat.indices for thermostat in temperature]...)
        if length(unique(indices)) != length(atoms.masses)
            throw(DomainError(temperature, "Every atom must have a TemperatureSetting applied to it."))
        end
        if length(indices) != length(unique(indices))
            throw(DomainError(temperature, "Atoms can only have one thermostat applied to them."))
        end
    end

    RingPolymerSimulation(temperature, cell, atoms, cache, DiabaticMDEF(atoms.masses, ndofs(model)), n_beads, solver)
end

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:Union{DiabaticMDEF,Classical},<:NQCCalculators.RingPolymer_QuantumModel_Cache}, t)
    
    NQCCalculators.update_cache!(sim.cache, r)
    adiabatic_derivative = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
    eigen = NQCCalculators.get_eigen(sim.cache, r)
    μ = NQCModels.fermilevel(sim.cache.model)
    β = 1 / get_temperature(sim, t)

    fill!(dv, zero(eltype(dv)))
    @views for i in axes(dv, 3) # bead
        NQCModels.state_independent_derivative!(sim.cache.model, dv[:,:,i], r[:,:,i])
        for j in axes(dv, 2) # atom
            for k in axes(dv, 1) # dof
                for n in eachindex(eigen[i].values) # state
                    ϵ = eigen[i].values[n]
                    f = fermi(ϵ, μ, β)
                    ∂f∂ϵ = ∂fermi(ϵ, μ, β)
                    ∂ϵ = adiabatic_derivative[k,j,i][n,n]
                    dv[k,j,i] += ∂ϵ*f + ∂f∂ϵ*∂ϵ*ϵ
                end
            end
        end
    end
    LinearAlgebra.lmul!(-1, dv)
    DynamicsUtils.divide_by_mass!(dv, masses(sim))

    return nothing
end

function evaluate_friction!(Λ::AbstractArray{T,3}, sim::RingPolymerSimulation{<:DiabaticMDEF}, r::AbstractArray{T,3}, t::Real) where {T}
    β = 1 / get_temperature(sim, t)
    μ = NQCModels.fermilevel(sim.cache.model)
    @views for i in NQCCalculators.beads(sim.cache)
        if sim.method.friction_method isa WideBandExact
            potential = NQCCalculators.get_potential(sim.cache, r)
            derivative = NQCCalculators.get_derivative(sim.cache, r)
            fill_friction_tensor!(Λ[:,:,i], sim.method.friction_method, potential[i], derivative[:,:,i], r[:,:,i], μ)
        else
            ∂H = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
            eigen = NQCCalculators.get_eigen(sim.cache, r)
            fill_friction_tensor!(Λ[:,:,i], sim.method.friction_method, ∂H[:,:,i], eigen[i], r[:,:,i], μ)
        end
    end
    return Λ
end

function DynamicsUtils.classical_potential_energy(
    sim::RingPolymerSimulation{<:Union{Classical,DiabaticMDEF},<:NQCCalculators.RingPolymer_QuantumModel_Cache},
    r::AbstractArray{T,3}
) where {T}
    eigen = NQCCalculators.get_eigen(sim.cache, r)
    μ = NQCModels.fermilevel(sim.cache.model)
    β = 1 / get_temperature(sim)

    potential = zero(T)
    @views for i in axes(r,3)
        potential += NQCModels.state_independent_potential(sim.cache.model, r[:,:,i])
        for ϵ in eigen[i].values
            potential += ϵ * fermi(ϵ, μ, β)
        end
    end
    return potential
end
