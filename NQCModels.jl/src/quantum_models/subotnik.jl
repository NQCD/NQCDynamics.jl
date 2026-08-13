"""
    MiaoSubotnik{T<:AbstractFloat} <: QuantumModel

Double well model with parameters matching those of Miao and Subotnik in the reference.
This model should be paired with the `AndersonHolstein` model to couple to the bath of metallic states.

# References

- J. Chem. Phys. 150, 041711 (2019)
"""
struct MiaoSubotnik{T<:AbstractFloat} <: QuantumModel
    Γ::T
    m::T
    ω::T
    g::T
    ΔG::T
    coupling::T
end

function MiaoSubotnik(;
    Γ,
    m = 2e3,
    ω = 2e-4,
    g = 20.6097,
    ΔG = -3.8e-3,
    coupling = sqrt(Γ/2π)
    )
    return MiaoSubotnik(Γ, m, ω, g, ΔG, coupling)
end

NQCModels.ndofs(::MiaoSubotnik) = 1
NQCModels.nstates(::MiaoSubotnik) = 2

function NQCModels.potential(model::MiaoSubotnik, R::AbstractMatrix)
    (;m, ω, g, ΔG, coupling) = model
    x = R[1]

    V11 = (m*ω^2*x^2)/2
    V22 = (m*ω^2*(x-g)^2)/2 + ΔG
    V12 = coupling

    return Hermitian([V11 V12; V12 V22])
end

function NQCModels.potential!(model::MiaoSubotnik, V::Hermitian, R::AbstractMatrix)
    (;m, ω, g, ΔG, coupling) = model
    x = R[1]

    V11 = (m*ω^2*x^2)/2
    V22 = (m*ω^2*(x-g)^2)/2 + ΔG
    V12 = coupling

    V.data .= Hermitian([V11 V12; V12 V22])
end

function NQCModels.derivative(model::MiaoSubotnik, R::AbstractMatrix)
    (;m, ω, g) = model
    x = R[1]

    D11 = m*ω^2*x
    D22 = m*ω^2*(x-g)
    
    return Hermitian([D11 0; 0 D22])
end

function NQCModels.derivative!(model::MiaoSubotnik, D::Hermitian, R::AbstractMatrix)
    (;m, ω, g) = model
    x = R[1]
    
    D11 = m*ω^2*x
    D22 = m*ω^2*(x-g)
    
    D.data .= [D11 0; 0 D22]
end
