"""
    MiaoSubotnik{T<:AbstractFloat} <: QuantumModel{Diabatic}

Double well model with parameters matching those of Miao and Subotnik in the reference.
This model should be paired with the `AndersonHolstein` model to couple to the bath of metallic states.

# References

- J. Chem. Phys. 150, 041711 (2019)
"""
struct MiaoSubotnik{T<:AbstractFloat} <: QuantumModel{Diabatic}
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

# GPU extensions: ext/NQCModelsCUDAExt.jl:254, ext/NQCModelsMetalExt.jl:254
function NQCModels.potential!(model::MiaoSubotnik, V::Hermitian, R::AbstractMatrix)
    (;m, ω, g, ΔG, coupling) = model

    V11 = sum(x -> (m*ω^2*x^2)/2, R)
    V22 = sum(x -> (m*ω^2*(x-g)^2)/2 + ΔG, R)
    V12 = coupling * length(R)

    V.data .= [V11 V12; V12 V22]
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:265, ext/NQCModelsMetalExt.jl:265
function NQCModels.derivative!(model::MiaoSubotnik, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    (;m, ω, g) = model

    for dof in axes(R, 1)
        for atom in axes(R, 2)
            x = R[dof,atom]
            D11 = m*ω^2*x
            D22 = m*ω^2*(x-g)

            D[dof,atom].data .= [D11 0; 0 D22]
        end
    end
end
