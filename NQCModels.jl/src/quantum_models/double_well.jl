
"""
    DoubleWell(mass=1, ω=1, γ=1, Δ=1)

Two state double well, also called the one-dimensional spin-boson model.
See: [J. Chem. Phys. 150, 244102 (2019)](https://doi.org/10.1063/1.5096276)
"""
Parameters.@with_kw struct DoubleWell{M,W,Y,D} <: QuantumModel
    mass::M = 1
    ω::W = 1
    γ::Y = 1
    Δ::D = 1
end

NQCModels.ndofs(::DoubleWell) = 1
NQCModels.nstates(::DoubleWell) = 2

function NQCModels.potential(model::DoubleWell, R::AbstractMatrix)
    r = R[1]

    V₀ = 0.5 * model.mass * model.ω^2 * r^2
    v = sqrt(2)*model.γ*r
    V11 = V₀ + v
    V22 = V₀ - v
    V12 = model.Δ/2

    return Hermitian([V11 V12; V12 V22])
end

function NQCModels.potential!(model::DoubleWell, V::Hermitian, R::AbstractMatrix)
    r = R[1]

    V₀ = 0.5 * model.mass * model.ω^2 * r^2
    v = sqrt(2)*model.γ*r
    V11 = V₀ + v
    V22 = V₀ - v
    V12 = model.Δ/2

    V.data .= [V11 V12; V12 V22]
end

# In-place derivative already defined in QuantumModels.jl

function NQCModels.derivative!(model::DoubleWell, D::Hermitian, R::AbstractMatrix)
    r = R[1]

    D₀ = model.mass * model.ω^2 * r
    v = sqrt(2)*model.γ
    D11 = D₀ + v
    D22 = D₀ - v
    D.data .= [D11 0; 0 D22]
end
