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

    @views for i in Calculators.beads(sim.calculator)
        NQCModels.state_independent_derivative!(sim.calculator.model, dv[:,:,i], r[:,:,i])
        current_bead_eigen = eigen[i]
        for j in NQCModels.mobileatoms(sim)
            for k in NQCModels.dofs(sim)
                for n in NQCModels.eachstate(sim)
                    f = fermi(current_bead_eigen.values[n], μ, β)
                    dv[k,j,i] += adiabatic_derivative[k,j,i][n,n] * f
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
            fill_friction_tensor!(Λ[:,:,i], sim.method.friction_method, potential[i], derivative[:,:,i], r[:,:,i], μ, β)
        else
            ∂H = Calculators.get_adiabatic_derivative(sim.calculator, r)
            eigen = Calculators.get_eigen(sim.calculator, r)
            fill_friction_tensor!(Λ[:,:,i], sim.method.friction_method, ∂H[:,:,i], eigen[i], r[:,:,i], μ, β)
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
    for (i, r_bead) in enumerate(eachbead(r))
        potential += NQCModels.state_independent_potential(sim.calculator.model, r_bead)
        for ϵ in eigen[i].values
            potential += ϵ * fermi(ϵ, μ, β)
        end
    end
    return potential
end
