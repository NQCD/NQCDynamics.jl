
struct Scattering1D <: QuantumFrictionModel
    n_states::UInt
    N::UInt
    a::Float64
    D::Float64
    α::Float64
    β::Float64
    B::Float64
    σ::Float64
end

NQCModels.ndofs(::Scattering1D) = 1
NQCModels.nstates(model::Scattering1D) = model.n_states

function Scattering1D(;N=10, a=1, D=1, α=0, β=1, B=1, σ=0.6u"eV")
    Scattering1D(N+1, N, a, D, α, β, B, austrip(σ))
end

function NQCModels.potential(model::Scattering1D, R::AbstractMatrix)
    r = R[1]

    V₀ = model.D*(exp(-2model.a*r) - 2*exp(-model.a*r))
    γᵣ = model.B*exp(-model.a*r^2)

    V = zeros((model.n_states, model.n_states))

    V[diagind(V,0)] .= model.α + V₀ #metal state energies
    V[diagind(V,1)] .= model.β #metal-metal couplings between neighbouring states
    
    V[end-1, end] = 0.0 #periodic coupling between first and last metal state
    V[2,end] = model.β 

    V[1,1] = V₀ # Occupied molecule state
    V[1,2] = γᵣ # molecule metal coupling

    return V
end

function NQCModels.potential!(model::Scattering1D, V::Hermitian, R::AbstractMatrix)
    r = R[1]

    V₀ = model.D*(exp(-2model.a*r) - 2*exp(-model.a*r))
    γᵣ = model.B*exp(-model.a*r^2)

    V.data[diagind(V,0)] .= model.α + V₀
    V.data[diagind(V,1)] .= model.β

    V.data[end-1, end] = 0.0
    V.data[2,end] = model.β

    V.data[1,1] = V₀ # Occupied molecule state
    V.data[1,2] = γᵣ # molecule metal coupling
end

function NQCModels.derivative!(model::Scattering1D, D::Hermitian, R::AbstractMatrix)
    r = R[1]
    
    D₀ = 2*model.D*model.a*(exp(-model.a*r)-exp(-2*model.a*r))
    dγᵣ = -2model.a*r*model.B*exp(-model.a*r^2)

    D.data[1,2] = dγᵣ
    D.data[diagind(D,0)] .= D₀
end
