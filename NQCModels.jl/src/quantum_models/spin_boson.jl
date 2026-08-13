using ..ClassicalModels: ClassicalModels

abstract type SpectralDensity end

"Discretize a given spectral density for N oscillators. Returns frequencies and couplings."
function discretize(J::SpectralDensity, N::Integer)
    ωⱼ = ω.(Ref(J), N, 1:N)
    cⱼ = c.(Ref(J), N, ωⱼ)
    return ωⱼ, cⱼ
end

"""
    OhmicSpectralDensity{T} <: SpectralDensity

Ohmic density as detailed in:
[Xin He, Jian Liu, J. Chem. Phys. 151, 024105 (2019)](http://aip.scitation.org/doi/10.1063/1.5108736)
"""
struct OhmicSpectralDensity{T} <: SpectralDensity
    ωᶜ::T
    α::T
end

(J::OhmicSpectralDensity)(ω) = π/2 * J.α * ω * exp(-ω / J.ωᶜ)
ω(J::OhmicSpectralDensity, N, j) = -J.ωᶜ * log(1 - j/(N+1))
c(J::OhmicSpectralDensity, N, ωⱼ) = sqrt(J.α * J.ωᶜ/(N+1)) * ωⱼ

"""
    DebyeSpectralDensity{T} <: SpectralDensity

Debye density as detailed in:
[Xin He, Jian Liu, J. Chem. Phys. 151, 024105 (2019)](http://aip.scitation.org/doi/10.1063/1.5108736)
"""
struct DebyeSpectralDensity{T} <: SpectralDensity
    ωᶜ::T
    λ::T
end

(J::DebyeSpectralDensity)(ω) = 2J.λ * J.ωᶜ * ω / (J.ωᶜ^2 + ω^2)
ω(J::DebyeSpectralDensity, N, j) = J.ωᶜ * tan(π/2 * (1 - j/(N+1)))
c(J::DebyeSpectralDensity, N, ωⱼ) = sqrt(2J.λ/(N+1)) * ωⱼ

"""
    AltDebyeSpectralDensity{T} <: SpectralDensity

Standard Debye spectral density but uses an alternative discretization scheme that requires
a cutoff parameter `ωᵐ`.

# References
[Najeh Rekik, Chang-Yu Hsieh, Holly Freedman, Gabriel Hanna, J. Chem. Phys. 138, 144106 (2013)](https://doi.org/10.1063/1.4799272)
"""
struct AltDebyeSpectralDensity{T} <: SpectralDensity
    ωᶜ::T
    λ::T
    ωᵐ::T
end

(J::AltDebyeSpectralDensity)(ω) = 2J.λ * J.ωᶜ * ω / (J.ωᶜ^2 + ω^2)
ω(J::AltDebyeSpectralDensity, N, j) = tan(j * atan(J.ωᵐ / J.ωᶜ) / N) * J.ωᶜ
c(J::AltDebyeSpectralDensity, N, ωⱼ) = sqrt(4J.λ * atan(J.ωᵐ / J.ωᶜ) / (π * N)) * ωⱼ

"""
    SpinBoson(density::SpectralDensity, N::Integer, ϵ, Δ)

Spin boson model with `N` bosons with given spectral density.

# References
[Xin He, Jian Liu, J. Chem. Phys. 151, 024105 (2019)](http://aip.scitation.org/doi/10.1063/1.5108736)
"""
struct SpinBoson{T} <: QuantumModel
    ϵ::T
    Δ::T
    ωⱼ::Vector{T}
    cⱼ::Vector{T}
end

NQCModels.nstates(::SpinBoson) = 2
NQCModels.ndofs(model::SpinBoson) = 1

function SpinBoson(density::SpectralDensity, N::Integer, ϵ, Δ)
    ωⱼ, cⱼ = discretize(density, N)
    SpinBoson(ϵ, Δ, ωⱼ, cⱼ)
end

function NQCModels.potential(model::SpinBoson, r::AbstractMatrix)

    Parameters.@unpack ωⱼ, cⱼ, ϵ, Δ = model

    temp = zeros(length(ωⱼ))'
    @. temp = ωⱼ'^2 * r^2 / 2

    v0 = sum(temp)

    V11 = v0 + ϵ
    V11 += sum(cⱼ' .* r)

    V22 = v0 - ϵ
    V22 -= sum(cⱼ' .* r)

    V12 = Δ

    return Hermitian([V11 V12; V12 V22])
end

function NQCModels.potential!(model::SpinBoson, V::Hermitian, r::AbstractMatrix)

    Parameters.@unpack ωⱼ, cⱼ, ϵ, Δ = model

    temp = zeros(length(ωⱼ))'
    @. temp = ωⱼ'^2 * r^2 / 2

    v0 = sum(temp)

    V11 = v0 + ϵ

    V11 += sum(cⱼ' .* r)
    
    V22 = v0 - ϵ
    V22 -= sum(cⱼ' .* r)

    V12 = Δ

    V.data .= [V11 V12; V12 V22]
end

function NQCModels.derivative(model::SpinBoson, R::AbstractMatrix)
    D = [Hermitian(zero(matrix_template(model, eltype(R)))) for _=1:size(R, 1), _=1:size(R, 2)]
    NQCModels.derivative!(model, D, R)
    return D
end

function NQCModels.derivative!(model::SpinBoson, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)

    Parameters.@unpack ωⱼ, cⱼ = model

    d0 = zero(r)
    @. d0 = ωⱼ'^2 * r

    for i in axes(r, 2)
        D[i].data .= [d0[1,i]+cⱼ[i] 0; 0 d0[1,i]-cⱼ[i]]
    end
end

"""
    BosonBath(density::SpectralDensity, N::Integer)

Bosonic bath with given spectral density.

Useful for sampling the bath uncoupled from the spin for spin-boson dynamics.
"""
struct BosonBath{T} <: ClassicalModels.ClassicalModel
    ωⱼ::Vector{T}
end

NQCModels.ndofs(model::BosonBath) = 1

function BosonBath(density::SpectralDensity, N::Integer)
    ωⱼ, _ = discretize(density, N)
    BosonBath(ωⱼ)
end

function NQCModels.potential(model::BosonBath, r::AbstractMatrix)
    return sum(model.ωⱼ'.^2 .* r.^2 ./2)
end

function NQCModels.potential!(model::BosonBath, V, r::AbstractMatrix)
    V = sum(model.ωⱼ'.^2 .* r.^2 ./2)
end

function NQCModels.derivative!(model::BosonBath, D::AbstractMatrix, r::AbstractMatrix)
    @. D = model.ωⱼ'^2 * r
end
