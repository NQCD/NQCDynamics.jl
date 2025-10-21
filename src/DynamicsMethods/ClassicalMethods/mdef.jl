using UnitfulAtomic
using NQCBase
using NQCModels
using NQCCalculators
using NQCDynamics: get_temperature, masses, TemperatureSetting
using Optim: Optim
using LinearAlgebra

abstract type AbstractMDEF <: DynamicsMethods.Method end

#= function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, friction!,
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end =#

#= function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:AbstractMDEF})
    StochasticDiffEq.DynamicalSDEProblem(acceleration!, DynamicsUtils.velocity!, SciMlBase.DynamicalNoiseFunction(matrix_friction_update!),
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end =#

"""
```math
dr = v dt\\\\
dv = -\\Delta U/M dt - \\Gamma v dt + \\sigma \\sqrt{2\\Gamma} dW
```
``\\Gamma`` is the friction tensor with units of inverse time.
For thermal dynamics we set ``\\sigma = \\sqrt{kT / M}``,
where ``T`` is the electronic temperature.

This is integrated using the BAOAB algorithm where the friction "O" step is performed
in the tensor's eigenbasis. See `src/dynamics/mdef_baoab.jl` for details.
"""
struct MDEF{M} <: AbstractMDEF
    mass_scaling::M
end

MDEF(masses::AbstractVector, DoFs::Integer) = MDEF(get_mass_scale_matrix(masses, DoFs))


struct DiabaticMDEF{M<:AbstractMatrix} <: AbstractMDEF
    mass_scaling::M
end

function DiabaticMDEF(masses::AbstractVector, DoFs::Integer)
    DiabaticMDEF(get_mass_scale_matrix(masses, DoFs))
end

function get_mass_scale_matrix(masses::AbstractVector, DoFs::Integer)
    vect = repeat(sqrt.(masses), inner=DoFs)
    vect * vect'
end

function NQCDynamics.Simulation{MDEF}(atoms::Atoms, model::Model; kwargs...)
    NQCDynamics.Simulation(atoms, model, MDEF(atoms.masses, ndofs(model));  kwargs...)
end
function NQCDynamics.RingPolymerSimulation{MDEF}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    NQCDynamics.RingPolymerSimulation(atoms, model, MDEF(atoms.masses, ndofs(model)), n_beads; kwargs...)
end

function NQCDynamics.Simulation{DiabaticMDEF}(atoms::Atoms{T}, model::Model; friction_method::FrictionEvaluationMethod, 
    temperature=0u"K", cell::AbstractCell=InfiniteCell(), solver::Symbol=:exact, kwargs...
    ) where {T}

    cache = Create_Cache(model, length(atoms), T; friction_method=friction_method)

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

    Simulation(temperature, cell, atoms, cache, DiabaticMDEF(atoms.masses, ndofs(model)), solver)
end

# Specific DiabaticMDEF Simulation multiple dispatch for the AndersonHolstein model
function NQCDynamics.Simulation{DiabaticMDEF}(atoms::Atoms{T}, model::AndersonHolstein; friction_method::FrictionEvaluationMethod, 
    temperature=0u"K", cell::AbstractCell=InfiniteCell(), solver::Symbol=:exact, kwargs...
    ) where {T}

    cache = Create_Cache(model, length(atoms), T, friction_method) # Create_Cache has for variable arguments, `friction_method`, is no-longer a keyword argument

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

    Simulation(temperature, cell, atoms, cache, DiabaticMDEF(atoms.masses, ndofs(model)), solver)
end

function acceleration!(dv, v, r, sim::Simulation{<:Union{DiabaticMDEF,Classical},<:NQCCalculators.Abstract_QuantumModel_Cache}, t)
    NQCCalculators.update_cache!(sim.cache, r)
    adiabatic_derivative = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
    eigen = NQCCalculators.get_eigen(sim.cache, r)
    μ = NQCModels.fermilevel(sim.cache.model)
    β = 1 / get_temperature(sim, t)

    fill!(dv, zero(eltype(dv)))
    NQCModels.state_independent_derivative!(sim.cache.model, dv, r)
    for I in eachindex(dv)
        for i in eachindex(eigen.values)
            ϵ = eigen.values[i]
            f = NQCCalculators.fermi(ϵ, μ, β)
            ∂f∂ϵ = NQCCalculators.∂fermi(ϵ, μ, β)
            ∂ϵ = adiabatic_derivative[I][i,i]
            dv[I] += ∂ϵ*f + ∂f∂ϵ*∂ϵ*ϵ
        end
    end
    LinearAlgebra.lmul!(-1, dv)
    DynamicsUtils.divide_by_mass!(dv, masses(sim))

    return nothing
end

"""
    friction!(g, r, sim, t)

Evaluates friction tensor
"""
function friction!(g, r, sim::AbstractSimulation{<:AbstractMDEF}, t)
    friction = NQCCalculators.get_friction(sim.cache, r)
    g .= friction ./ sim.method.mass_scaling
    g
end

function friction!(g, r, sim::AbstractSimulation{<:DiabaticMDEF}, t)
    friction = NQCCalculators.get_friction(sim.cache, r)
    g .= friction ./ sim.method.mass_scaling
end


function determine_fermi_level(nelectrons, β, eigenvalues)

    count_electrons(μ) = sum(NQCCalculators.fermi(ϵ, μ, β) for ϵ in eigenvalues)
    optim_func(μ) = (count_electrons(μ) - nelectrons)^2

    optim = Optim.optimize(optim_func, eigenvalues[begin], eigenvalues[end])
    μ = Optim.minimizer(optim)
    filled_electrons = count_electrons(μ)
    if !isapprox(filled_electrons, nelectrons)
        throw(error(
            "Unable to determine the fermi level. \
            Got $μ with $filled_electrons electrons but there should be $nelectrons electrons. \
            Try increasing the temperature."
        ))
    end

    return μ
end

function DynamicsUtils.classical_potential_energy(
    sim::Simulation{<:Union{DiabaticMDEF,Classical},<:NQCCalculators.Abstract_QuantumModel_Cache},
    r::AbstractMatrix
)
    eigen = NQCCalculators.get_eigen(sim.cache, r)
    μ = NQCModels.fermilevel(sim.cache.model)
    β = 1 / get_temperature(sim)

    potential = NQCModels.state_independent_potential(sim.cache.model, r)
    for ϵ in eigen.values
        potential += ϵ * NQCCalculators.fermi(ϵ, μ, β)
    end
    return potential
end

function matrix_friction_update!(dutmp, utmp, integrator, cache)
    int = integrator

    DynamicsMethods.ClassicalMethods.friction!(cache.gtmp, int.utmp, int.p, int.t+int.dt*cache.half)
    DynamicsMethods.IntegrationAlgorithms.step_O!(cache, int)
end