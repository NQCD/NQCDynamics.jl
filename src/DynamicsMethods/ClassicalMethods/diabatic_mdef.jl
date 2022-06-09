using NQCModels: NQCModels
using NQCDynamics: get_temperature, masses
using Optim: Optim
using QuadGK: QuadGK
using LinearAlgebra: LinearAlgebra, diagind

abstract type FrictionEvaluationMethod end

struct DiabaticMDEF{M<:AbstractMatrix,F<:FrictionEvaluationMethod} <: AbstractMDEF
    mass_scaling::M
    friction_method::F
end

function DiabaticMDEF(masses::AbstractVector, DoFs::Integer, friction_method)
    DiabaticMDEF(get_mass_scale_matrix(masses, DoFs), friction_method)
end

function NQCDynamics.Simulation{DiabaticMDEF}(atoms::Atoms, model::Model;
    friction_method, kwargs...
)
    NQCDynamics.Simulation(atoms, model,
        DiabaticMDEF(atoms.masses, ndofs(model), friction_method);
        kwargs...
    )
end

function acceleration!(dv, v, r, sim::Simulation{<:DiabaticMDEF}, t)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    eigen = Calculators.get_eigen(sim.calculator, r)
    μ = NQCModels.fermilevel(sim.calculator.model)
    β = 1 / get_temperature(sim, t)

    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)
    for I in eachindex(dv)
        for i in eachindex(eigen.values)
            f = fermi(eigen.values[i], μ, β)
            dv[I] -= adiabatic_derivative[I][i,i] * f
        end
    end
    DynamicsUtils.divide_by_mass!(dv, masses(sim))

    return nothing
end

function friction!(g, r, sim::Simulation{<:DiabaticMDEF}, t)
    evaluate_friction!(g, sim, r, t)
    g ./= sim.method.mass_scaling
end

fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ-μ)))
function ∂fermi(ϵ, μ, β)
    ∂f = -β * exp(β*(ϵ-μ)) / (1 + exp(β*(ϵ-μ)))^2
    return isnan(∂f) ? zero(ϵ) : ∂f
end

gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))
gauss(x, friction_method::FrictionEvaluationMethod) = gauss(x, friction_method.σ)

function evaluate_friction!(Λ::AbstractMatrix, sim::Simulation{<:DiabaticMDEF}, r::AbstractMatrix, t::Real)
    β = 1/get_temperature(sim, t)
    μ = NQCModels.fermilevel(sim.calculator.model)
    fill_friction_tensor!(Λ, sim.method.friction_method, sim.calculator, r, μ, β)
    return Λ
end

function fill_friction_tensor!(Λ, friction_method::FrictionEvaluationMethod, calculator, r, μ, β)
    ∂H = Calculators.get_adiabatic_derivative(calculator, r)
    eigen = Calculators.get_eigen(calculator, r)
    for I in eachindex(r)
        for J in eachindex(r)
            Λ[J,I] = friction_method(∂H[J], ∂H[I], eigen.values, μ, β)
        end
    end
end

struct GaussianBroadening{T} <: FrictionEvaluationMethod 
    σ::T
end
function (friction_method::GaussianBroadening)(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        for m in eachindex(eigenvalues)
            ϵₙ = eigenvalues[n]
            ϵₘ = eigenvalues[m]
            Δϵ = ϵₙ - ϵₘ
            out += -π * ∂Hᵢ[n,m] * ∂Hⱼ[m,n] * gauss(Δϵ, friction_method) * ∂fermi(ϵₙ, μ, β)
        end
    end
    return out
end

struct OffDiagonalGaussianBroadening{T} <: FrictionEvaluationMethod
    σ::T
end
function (friction_method::OffDiagonalGaussianBroadening)(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        for m=n+1:length(eigenvalues)
            ϵₙ = eigenvalues[n]
            ϵₘ = eigenvalues[m]
            Δϵ = ϵₙ - ϵₘ

            fₙ = fermi(ϵₙ, μ, β)
            fₘ = fermi(ϵₘ, μ, β)
            Δf = (fₘ - fₙ)

            out += 2π * ∂Hᵢ[n,m] * ∂Hⱼ[m,n] * gauss(Δϵ, friction_method) * Δf / Δϵ
        end
    end
    return out
end

struct DirectQuadrature{T} <: FrictionEvaluationMethod
    ρ::T    
end
function (friction_method::DirectQuadrature)(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        ϵₙ = eigenvalues[n]
        out += -π * ∂Hᵢ[n,n] * ∂Hⱼ[n,n] * friction_method.ρ * ∂fermi(ϵₙ, μ, β)
    end
    return out
end

struct WideBandExact{T} <: FrictionEvaluationMethod
    ρ::T
end
function (friction_method::WideBandExact)(potential, ∂potentialᵢ, ∂potentialⱼ, μ, β)
    h = potential[1,1]
    ∂hᵢ = ∂potentialᵢ[1,1]
    ∂hⱼ = ∂potentialⱼ[1,1]

    ρ = friction_method.ρ
    Γ = 2π * potential[2,1]^2 * ρ
    ∂Γ∂potential = sqrt(8π * ρ * Γ)
    ∂Γᵢ = ∂Γ∂potential * ∂potentialᵢ[2,1]
    ∂Γⱼ = ∂Γ∂potential * ∂potentialⱼ[2,1]

    A(ϵ) = 1/π * Γ/2 / ((ϵ-h)^2 + (Γ/2)^2)
    kernel(ϵ) = -π * (∂hᵢ + (ϵ-h)*∂Γᵢ/Γ) * (∂hⱼ + (ϵ-h)*∂Γⱼ/Γ) * A(ϵ)^2 * ∂fermi(ϵ, μ, β)
    diagonal = @view potential[diagind(potential)]
    integral, _ = QuadGK.quadgk(kernel, extrema(diagonal)...)
    return integral
end

function fill_friction_tensor!(Λ, friction_method::WideBandExact, calculator, r, μ, β)
    potential = Calculators.get_potential(calculator, r)
    derivative = Calculators.get_derivative(calculator, r)
    for I in eachindex(r)
        for J in eachindex(r)
            Λ[J,I] = friction_method(potential, derivative[J], derivative[I], μ, β)
        end
    end
end

function determine_fermi_level(nelectrons, β, eigenvalues)

    count_electrons(μ) = sum(fermi(ϵ, μ, β) for ϵ in eigenvalues)
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

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:DiabaticMDEF}, r::AbstractMatrix)
    eigen = Calculators.get_eigen(sim.calculator, r)
    μ = NQCModels.fermilevel(sim.calculator.model)
    β = 1 / get_temperature(sim)

    potential = NQCModels.state_independent_potential(sim.calculator.model, r)
    for ϵ in eigen.values
        potential += ϵ * fermi(ϵ, μ, β)
    end
    return potential
end
