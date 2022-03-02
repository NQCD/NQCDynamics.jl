using NQCModels: NQCModels
using NQCDynamics: get_temperature
using Optim: Optim

struct DiabaticMDEF{T,M} <: AbstractMDEF
    mass_scaling::M
    σ::T
    friction_type::Symbol
end

function DiabaticMDEF(masses::AbstractVector, DoFs::Integer, σ, friction_type)
    DiabaticMDEF(get_mass_scale_matrix(masses, DoFs), austrip(σ), friction_type)
end

function NQCDynamics.Simulation{DiabaticMDEF}(atoms::Atoms, model::Model;
    σ=1.0, friction_type=:GB,
    kwargs...
)
    NQCDynamics.Simulation(atoms, model,
        DiabaticMDEF(atoms.masses, ndofs(model), σ, friction_type);
        kwargs...
    )
end

function acceleration!(dv, v, r, sim::Simulation{<:DiabaticMDEF}, t)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    eigen = Calculators.get_eigen(sim.calculator, r)
    μ = NQCModels.fermilevel(sim.calculator.model)
    β = 1 / get_temperature(sim.calculator.model, t)

    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)
    for I in eachindex(dv)
        for i in eachindex(eigen.values)
            f = fermi(eigen.values[i], μ, β)
            dv[I] -= adiabatic_derivative[I][i,i] * f
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masse)

    return nothing
end

function friction!(g, r, sim::Simulation{<:DiabaticMDEF}, t)
    # Calculators.evaluate_friction!(sim.calculator, r)
    # g .= sim.calculator.friction ./ sim.method.mass_scaling
end

fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ-μ)))
function ∂fermi(ϵ, μ, β)
    ∂f = -β * exp(β*(ϵ-μ)) / (1 + exp(β*(ϵ-μ)))^2
    return isnan(∂f) ? zero(ϵ) : ∂f
end
gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))

function evaluate_friction!(Λ::AbstractMatrix, sim::Simulation{<:DiabaticMDEF}, r::AbstractMatrix, t::Real)
    ∂H = Calculators.get_adiabatic_derivative(sim.calculator, r)
    eigen = Calculators.get_eigen(sim.calculator, r)
    σ = sim.method.σ
    β = 1/get_temperature(sim, t)
    μ = NQCModels.fermilevel(sim.calculator.model)

    ρ = sim.calculator.model.ρ
    potential = Calculators.get_potential(sim.calculator, r)
    derivative = Calculators.get_derivative(sim.calculator, r)

    for I in eachindex(r)
        for J in eachindex(r)
            if sim.method.friction_type === :GB
                Λ[I,J] = friction_gaussian_broadening(∂H[I], ∂H[J], eigen.values, μ, β, σ)
            elseif sim.method.friction_type === :ONGB
                Λ[I,J] = friction_off_diagonal_gaussian_broadening(∂H[I], ∂H[J], eigen.values, μ, β, σ)
            elseif sim.method.friction_type === :DQ
                Λ[I,J] = friction_direct_quadrature(∂H[I], ∂H[J], eigen.values, μ, β, ρ)
            elseif sim.method.friction_type === :WB
                Λ[I,J] = friction_wideband(potential, derivative[I], eigen.values, μ, β, ρ)
            else
                throw(ArgumentError("Friction type $(sim.method.friction_type) not recognised."))
            end
        end
    end
end

function friction_gaussian_broadening(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β, σ)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        for m in eachindex(eigenvalues)
            ϵₙ = eigenvalues[n]
            ϵₘ = eigenvalues[m]
            Δϵ = ϵₙ - ϵₘ
            out += -π * ∂Hᵢ[n,m] * ∂Hⱼ[m,n] * gauss(Δϵ, σ) * ∂fermi(ϵₙ, μ, β)
        end
    end
    return out
end

function friction_off_diagonal_gaussian_broadening(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β, σ)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        for m=n+1:length(eigenvalues)
            ϵₙ = eigenvalues[n]
            ϵₘ = eigenvalues[m]
            Δϵ = ϵₙ - ϵₘ

            fₙ = fermi(ϵₙ, μ, β)
            fₘ = fermi(ϵₘ, μ, β)
            Δf = (fₘ - fₙ)

            out += 2π * ∂Hᵢ[n,m] * ∂Hⱼ[m,n] * gauss(Δϵ, σ) * Δf / Δϵ
        end
    end
    return out
end

function friction_direct_quadrature(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β, ρ)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        ϵₙ = eigenvalues[n]
        out += -π * ∂Hᵢ[n,n] * ∂Hⱼ[n,n] * ρ * ∂fermi(ϵₙ, μ, β)
    end
    return out
end

using QuadGK: QuadGK

function friction_wideband(potential, derivative, eigenvalues, μ, β, ρ)
    h = potential[1,1]
    ∂h = derivative[1,1]
    Γ = 2π * potential[2,1]^2 * ρ
    ∂Γ∂potential = 4π * potential[2,1] * ρ
    ∂Γ = ∂Γ∂potential * derivative[2,1]

    A(ϵ) = 1/π * Γ/2 / ((ϵ-h)^2 + (Γ/2)^2)
    kernel(ϵ) = -π * (∂h + (ϵ-h)*∂Γ/Γ)^2 * A(ϵ)^2 * ∂fermi(ϵ, μ, β)
    integral, _ = QuadGK.quadgk(kernel, eigenvalues[begin], eigenvalues[end])
    return integral
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

