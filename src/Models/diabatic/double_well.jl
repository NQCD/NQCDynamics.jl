export DoubleWell

"""
    DoubleWell(mass=1, ω=1, γ=1, Δ=1)

Two state double well, also called the one-dimensional spin-boson model.
See: [J. Chem. Phys. 150, 244102 (2019)](https://doi.org/10.1063/1.5096276)
"""
@with_kw struct DoubleWell{M,W,Y,D} <: DiabaticModel
    n_states::UInt8 = 2
    mass::M = 1
    ω::W = 1
    γ::Y = 1
    Δ::D = 1
end

function potential!(model::DoubleWell, V::Hermitian, R::AbstractMatrix)

    V0(R) = 0.5 * model.mass * model.ω^2 * R^2

    V₀ = V0(R[1])
    v = sqrt(2)*model.γ*R[1] # Access only R[1] as this is a 1D model
    V[1,1] = V₀ + v
    V[2,2] = V₀ - v
    V.data[1,2] = model.Δ/2 # Sets both off-diagonals as V is Hermitian
end

function derivative!(model::DoubleWell, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)

    D0(R) = model.mass * model.ω^2 * R

    D₀ = D0(R[1])
    v = sqrt(2)*model.γ
    D[1][1,1] = D₀ + v
    D[1][2,2] = D₀ - v
end
