export Metal
using Unitful
using UnitfulAtomic

struct Metal{T} <: DiabaticFrictionModel
    n_states::Int
    N::Int
    mass::T
    ω::T
    γ::T
    Δ::T
    β::T
end

Metal(N; mass=1, ω=1, γ=1, Δ=1, β=0) = Metal{Float64}(N, N, mass, ω, γ, Δ, β)

# function Scattering1D(;N=10, a=1, D=1, α=0, β=1, B=1, σ=0.6u"eV")
#     Scattering1D(N+1, N, a, D, α, β, B, austrip(σ))
# end

function v1(model::Metal, V::Hermitian, R::AbstractMatrix)

    V0(R) = 0.5 * model.mass * model.ω^2 * R^2


    V₀ = V0(R[1])
    v = sqrt(2)*model.γ*R[1] # Access only R[1] as this is a 1D model
    # V[1,1] = V₀ + v
    # V[2,2] = V₀ - v
    # V.data[1,2] = model.Δ/2 # Sets both off-diagonals as V is Hermitian
    Γ = R[1]
    V[diagind(V)] .= V₀

    # V[1,1] += v
    # V[2,2] -= v
    # V.data[1,2] = model.Δ/2

    ϵ = range(0, 1, length=2model.N)
    # for i in 3:model.n_states
    #     V[i,i] += ϵ[i-2]
    #     V.data[1,i] = model.β  
    #     V.data[2,i] = 2model.β  
    # end
    for i in 1:2:2model.N
        V[i,i] += v + ϵ[i]
        V[i+1,i+1] -= v - ϵ[i]
        V.data[i,i+1] = model.Δ/2
    end

    for i in 3:2:2model.N
        V[i,i] -= V₀ + v
        V[i+1,i+1] -= V₀ - v
        V.data[1, i] = model.β
        V.data[2, i+1] = model.β
    end

end

function potential!(model::Metal, V::Hermitian, R::AbstractMatrix)

    v0(x) = model.mass * model.ω^2 * (x - model.γ/2)^2 / 2
    v1(x) = model.mass * model.ω^2 * (x + model.γ/2)^2 / 2

    V₀ = v0(R[1])
    V[diagind(V)] .= V₀

    V[2,2] = v1(R[1])

    ϵ = range(0, 2, length=model.n_states)
    for i=2:model.n_states
        V[i,i] += ϵ[i-1]
        V.data[1,i] = model.β
    end

end

# function derivative!(model::Metal, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
#     D0(R) = 2*model.D*model.a*(exp(-model.a*R)-exp(-2*model.a*R))
#     dγ(R) = -2model.a*R*model.B*exp(-model.a*R^2)

#     D₀ = D0(R[1])
#     dγᵣ = dγ(R[1])

#     D[1].data[1,2] = dγᵣ
#     for i=1:model.N+1
#         D[1][i,i] = D₀
#     end
# end
