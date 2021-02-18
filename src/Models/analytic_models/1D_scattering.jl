export Scattering1D
using Unitful
using UnitfulAtomic

struct Scattering1D <: DiabaticFrictionModel
    n_states::UInt
    N::UInt
    a::Float64
    D::Float64
    α::Float64
    β::Float64
    B::Float64
    σ::Float64
end

function Scattering1D(;N=10, a=1, D=1, α=0, β=1, B=1, σ=0.6u"eV")
    Scattering1D(N+1, N, a, D, α, β, B, austrip(σ))
end

function potential!(model::Scattering1D, V::Hermitian, R::AbstractMatrix)
    V0(R) = model.D*(exp(-2model.a*R) - 2*exp(-model.a*R))
    γ(R) = model.B*exp(-model.a*R^2)

    V₀ = V0(R[1])
    γᵣ = γ(R[1])

    V[1,1] = V₀ # Occupied molecule state
    V.data[1,2] = γᵣ # molecule metal coupling
    for i=2:model.N
        V[i,i] = model.α + V₀
        V.data[i,i+1] = model.β
    end
    V.data[2,end] = model.β
    V[end,end] = model.α + V₀
end

function derivative!(model::Scattering1D, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    D0(R) = 2*model.D*model.a*(exp(-model.a*R)-exp(-2*model.a*R))
    dγ(R) = -2model.a*R*model.B*exp(-model.a*R^2)

    D₀ = D0(R[1])
    dγᵣ = dγ(R[1])

    D[1].data[1,2] = dγᵣ
    for i=1:model.N+1
        D[1][i,i] = D₀
    end
end
