abstract type FrictionEvaluationMethod end

struct WideBandExact{T} <: FrictionEvaluationMethod
    ρ::T
    β::T
end

struct GaussianBroadening{T} <: FrictionEvaluationMethod 
    σ::T
    β::T
end

struct OffDiagonalGaussianBroadening{T} <: FrictionEvaluationMethod
    σ::T
    β::T
end

struct DirectQuadrature{T} <: FrictionEvaluationMethod
    ρ::T    
    β::T
end

#Definitions for fill_friction_tensor!()

"""
    function: fill_friction_tensor!(Λ, friction_method, ∂H, eigen, r, μ)

Λ is a friction tensor that is populated by the function using the friction evaluation method passed into the function.
∂H, eigen, r, μ are variables used by the friction method to evaluate Λ.
"""
function fill_friction_tensor!(Λ, friction_method::FrictionEvaluationMethod, eigen, ∂H, r, μ)
    for I in eachindex(r)
        for J in eachindex(r)
            Λ[J,I] = friction_method(∂H[J], ∂H[I], eigen.w, μ, friction_method.β)
        end
    end
end

"""
    function: fill_friction_tensor!(Λ, friction_method, potential, derivative, r, μ)

Λ is a friction tensor that is populated by the function using the wide-band approximation.
∂H, eigen, r, μ are variables used by the friction method to evaluate Λ.
"""
function fill_friction_tensor!(Λ, friction_method::WideBandExact, potential, derivative, r, μ)
    for J in eachindex(r)
        for I in eachindex(r)
            Λ[I,J] = friction_method(potential, derivative[I], derivative[J], μ, friction_method.β)
        end
    end
end

#Definitions of friction method functions

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
    integral, _ = QuadGK.quadgk(kernel, -Inf, Inf; rtol=1e-6)
    return integral
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

function (friction_method::DirectQuadrature)(∂Hᵢ, ∂Hⱼ, eigenvalues, μ, β)
    out = zero(eltype(eigenvalues))
    for n in eachindex(eigenvalues)
        ϵₙ = eigenvalues[n]
        out += -π * ∂Hᵢ[n,n] * ∂Hⱼ[n,n] * friction_method.ρ * ∂fermi(ϵₙ, μ, β)
    end
    return out
end

#Convenience functions

fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ-μ)))

function ∂fermi(ϵ, μ, β)
    ∂f = -β * exp(β*(ϵ-μ)) / (1 + exp(β*(ϵ-μ)))^2
    return isnan(∂f) ? zero(ϵ) : ∂f
end

gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))
gauss(x, friction_method::FrictionEvaluationMethod) = exp(-0.5 * x^2 / friction_method.σ^2) / (friction_method.σ*sqrt(2π))
