using NQCDynamics: nbeads
using RingPolymerArrays: eachbead

function NQCDynamics.RingPolymerSimulation{DiabaticMDEF}(atoms::Atoms, model::Model, n_beads::Integer;
    friction_method, kwargs...
)
    NQCDynamics.RingPolymerSimulation(atoms, model,
        DiabaticMDEF(atoms.masses, ndofs(model), friction_method),
        n_beads;
        kwargs...
    )
end

function acceleration!(dv, v, r, sim::RingPolymerSimulation{<:Union{DiabaticMDEF,Classical},<:Calculators.RingPolymerLargeDiabaticCalculator}, t)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    eigen = Calculators.get_eigen(sim.calculator, r)
    μ = NQCModels.fermilevel(sim.calculator.model)
    β = 1 / get_temperature(sim, t)

    fill!(dv, zero(eltype(dv)))
    @views for i in axes(dv, 3) # bead
        NQCModels.state_independent_derivative!(sim.calculator.model, dv[:,:,i], r[:,:,i])
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
    μ = NQCModels.fermilevel(sim.calculator.model)
    @views for i in Calculators.beads(sim.calculator)
        if sim.method.friction_method isa WideBandExact
            potential = Calculators.get_potential(sim.calculator, r)
            derivative = Calculators.get_derivative(sim.calculator, r)
            fill_friction_tensor!(Λ[:,:,i], sim.method.friction_method, potential[i], derivative[:,:,i], r[:,:,i], μ)
        else
            ∂H = Calculators.get_adiabatic_derivative(sim.calculator, r)
            eigen = Calculators.get_eigen(sim.calculator, r)
            fill_friction_tensor!(Λ[:,:,i], sim.method.friction_method, ∂H[:,:,i], eigen[i], r[:,:,i], μ)
        end
    end
    return Λ
end

function DynamicsUtils.classical_potential_energy(
    sim::RingPolymerSimulation{<:Union{Classical,DiabaticMDEF},<:Calculators.RingPolymerLargeDiabaticCalculator},
    r::AbstractArray{T,3}
) where {T}
    eigen = Calculators.get_eigen(sim.calculator, r)
    μ = NQCModels.fermilevel(sim.calculator.model)
    β = 1 / get_temperature(sim)

    potential = zero(T)
    @views for i in axes(r,3)
        potential += NQCModels.state_independent_potential(sim.calculator.model, r[:,:,i])
        for ϵ in eigen[i].values
            potential += ϵ * fermi(ϵ, μ, β)
        end
    end
    return potential
end
