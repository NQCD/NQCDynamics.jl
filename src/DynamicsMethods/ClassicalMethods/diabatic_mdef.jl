using NQCModels: NQCModels
using NQCDynamics: get_temperature

struct DiabaticMDEF{T,M} <: AbstractMDEF
    mass_scaling::M
    fermi_level::T
    σ::T
    friction_type::Symbol
end

function DiabaticMDEF(masses::AbstractVector, DoFs::Integer, fermi_level, σ, friction_type)
    DiabaticMDEF(get_mass_scale_matrix(masses, DoFs), austrip(fermi_level), austrip(σ), friction_type)
end

function NQCDynamics.Simulation{DiabaticMDEF}(atoms::Atoms, model::Model;
    fermi_level=0.0, σ=1.0, friction_type=:GB,
    kwargs...
)
    NQCDynamics.Simulation(atoms, model,
        DiabaticMDEF(atoms.masses, ndofs(model), fermi_level, σ, friction_type);
        kwargs...
    )
end

function acceleration!(dv, v, r, sim::Simulation{<:DiabaticMDEF}, t)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    eigen = Calculators.get_eigen(sim.calculator, r)

    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)
    for I in eachindex(dv)
        for i in eachindex(eigen.values)
            f = fermi(eigen.values[i], sim.method.fermi_level, 1/get_temperature(sim, t))
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
∂fermi(ϵ, μ, β) = -β * exp(β*(ϵ-μ)) / (1 + exp(β*(ϵ-μ)))^2
gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))

function evaluate_friction!(Λ::AbstractMatrix, sim::Simulation{<:DiabaticMDEF}, r::AbstractMatrix, t::Real)
    ∂H = Calculators.get_adiabatic_derivative(sim.calculator, r)
    eigen = Calculators.get_eigen(sim.calculator, r)
    μ = NQCModels.fermilevel(sim.calculator.model)
    σ = sim.method.σ
    β = 1/get_temperature(sim, t)

    fill!(Λ, zero(eltype(r)))
    for I in eachindex(r)
        for J in eachindex(r)
            if sim.method.friction_type === :GB
                Λ[I,J] = friction_gaussian_broadening(∂H[I], ∂H[J], eigen.values, μ, β, σ)
            elseif sim.method.friction_type === :ONGB
                Λ[I,J] = friction_off_diagonal_gaussian_broadening(∂H[I], ∂H[J], eigen.values, μ, β, σ)
            elseif sim.method.friction_type === :DQ
                ρ = 1 / (sim.calculator.model.bathstates[1] - sim.calculator.model.bathstates[2])
                Λ[I,J] = friction_direct_quadrature(∂H[I], ∂H[J], eigen.values, μ, β, ρ)
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
            # if n != m
                ϵₙ = eigenvalues[n]
                ϵₘ = eigenvalues[m]
                Δϵ = ϵₙ - ϵₘ
                out += -π * ∂Hᵢ[n,m] * ∂Hⱼ[m,n] * gauss(Δϵ, σ) * ∂fermi(ϵₙ, μ, β)
            # end
        end
    end
    return out
end

function friction_off_diagonal_gaussian_broadening(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β, σ)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        # for m in eachindex(eigenvalues)
        for m=n+1:length(eigenvalues)
            if n != m
                ϵₙ = eigenvalues[n]
                ϵₘ = eigenvalues[m]
                Δϵ = ϵₙ - ϵₘ

                fₙ = fermi(ϵₙ, μ, β)
                fₘ = fermi(ϵₘ, μ, β)
                Δf = (fₘ - fₙ)

                out += 2π * ∂Hᵢ[n,m] * ∂Hⱼ[m,n] * gauss(Δϵ, σ) * Δf / Δϵ
            end
        end
    end
    return out
end

function friction_direct_quadrature(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β, ρ)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        ϵₙ = eigenvalues[n]
        out += π * ∂Hᵢ[n,n] * ∂Hⱼ[n,n] * ρ * ∂fermi(ϵₙ, μ, β)
    end
    return out
end

