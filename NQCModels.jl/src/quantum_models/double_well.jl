
"""
    DoubleWell(mass=1, ω=1, γ=1, Δ=1)

Two state double well, also called the one-dimensional spin-boson model.
See: [J. Chem. Phys. 150, 244102 (2019)](https://doi.org/10.1063/1.5096276)
"""
Parameters.@with_kw struct DoubleWell{M,W,Y,D} <: QuantumModel{Diabatic}
    mass::M = [1.0]
    ω::W = [1.0]
    γ::Y = [1.0]
    Δ::D = [1.0]
end

NQCModels.ndofs(::DoubleWell) = 1
NQCModels.nstates(::DoubleWell) = 2

# GPU extensions: ext/NQCModelsCUDAExt.jl:180, ext/NQCModelsMetalExt.jl:180
function NQCModels.potential(model::DoubleWell, R::AbstractMatrix)
    V₀ = @. 0.5 * model.mass' * model.ω^2 * R^2
    v = @. sqrt(2)*model.γ*R
    V11 = @. V₀ + v
    V22 = @. V₀ - v
    V12 = @. model.Δ / 2 + 0 * R
    return Hermitian([sum(V11) sum(V12); sum(V12) sum(V22)])
end

function NQCModels.potential!(model::DoubleWell, V::Hermitian, R::AbstractMatrix)
    V.data .= NQCModels.potential(model, R)
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:192, ext/NQCModelsMetalExt.jl:192
function NQCModels.derivative!(model::DoubleWell, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    getparam(v, atom) = length(v) == 1 ? v[1] : v[atom]
    for dof in axes(r, 1)
        for atom in axes(r, 2)
            mass = getparam(model.mass, atom)
            ω = getparam(model.ω, atom)
            γ = getparam(model.γ, atom)

            D₀ = mass * ω^2 * r[1,atom]
            v = sqrt(2)*γ
            D11 = D₀ + v
            D22 = D₀ - v
            D[dof,atom].data .= [D11 0; 0 D22]
        end
    end
end
