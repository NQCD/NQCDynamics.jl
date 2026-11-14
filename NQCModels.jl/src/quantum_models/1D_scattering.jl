# Currently multi-particle capable, but all with same parameters.
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

# GPU extensions: ext/NQCModelsCUDAExt.jl:439, ext/NQCModelsMetalExt.jl:439
function NQCModels.potential!(model::Scattering1D, V::Hermitian, R::AbstractMatrix)
    M = length(R)
    V₀(r) = model.D*(exp(-2model.a*r) - 2*exp(-model.a*r))

    V.data .= 0
    V.data[diagind(V,0)[2:end]] .= sum(r -> model.α + V₀(r), R) # metal state energies
    V.data[diagind(V,1)] .= M * model.β # metal-metal couplings between neighbouring states

    V.data[end-1, end] = 0.0 # periodic coupling between first and last metal state
    V.data[2,end] = M * model.β

    V.data[1,1] = sum(V₀, R) # Occupied molecule state
    V.data[1,2] = sum(r -> model.B*exp(-model.a*r^2), R) # molecule metal coupling
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:459, ext/NQCModelsMetalExt.jl:459
function NQCModels.derivative!(model::Scattering1D, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    for atom in axes(R, 2) # Hardcoded to single DoF
        r = R[1, atom]
        D₀ = 2*model.D*model.a*(exp(-model.a*r)-exp(-2*model.a*r))
        dγᵣ = -2model.a*r*model.B*exp(-model.a*r^2)

        D[1,atom].data[1,2] = dγᵣ
        D[1,atom].data[diagind(D[1,atom],0)] .= D₀
    end
end
