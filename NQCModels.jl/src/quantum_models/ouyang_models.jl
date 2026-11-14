using LinearAlgebra: diagind

"""
    OuyangModelOne(A=0.01, B=1.6, Γ=1e-4, N=10, ΔE=1.6e-2, D=1.0)

Model #1 from [Ouyang and Subotnik](http://aip.scitation.org/doi/10.1063/1.4908032).
See also [Ouyang's thesis](https://repository.upenn.edu/edissertations/2508/).
"""
Parameters.@with_kw struct OuyangModelOne <: QuantumFrictionModel
    A::Float64 = 0.01
    B::Float64 = 1.6
    Γ::Float64 = 1e-4
    N::Int = 10
    ΔE::Float64 = 1.6e-2
    ρ::Float64 = N / ΔE
    C::Float64 = sqrt(Γ / 2π*ρ) / N
    D::Float64 = 1.0

    n_states::Int = N + 1
end

NQCModels.nstates(model::OuyangModelOne) = model.n_states
NQCModels.ndofs(::OuyangModelOne) = 1

# GPU extensions: ext/NQCModelsCUDAExt.jl:409, ext/NQCModelsMetalExt.jl:409
function NQCModels.potential!(model::OuyangModelOne, V::Hermitian, R::AbstractMatrix)
    Parameters.@unpack A, B, C, D, ΔE, N = model
    Δs = range(-ΔE/2, ΔE/2, length=N)
    M = length(R)

    V.data .= 0
    V.data[1,1] = sum(r -> A*tanh(B*r), R)
    V.data[1,2:N+1] .= sum(r -> C*exp(-D*r^2), R)
    V.data[diagind(N+1,N+1)[2:end]] .= sum(r -> -(A*tanh(B*r)), R) .+ M .* Δs
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:424, ext/NQCModelsMetalExt.jl:424
function NQCModels.derivative!(model::OuyangModelOne, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    (;A, B, C, N) = model
    Dparam = model.D

    for atom in axes(R, 2)
        r = R[1, atom]
        derivative = D[1,atom]
        derivative.data[1,1] = B*A*sech(B*r)^2
        derivative.data[1,2:N+1] .= -2Dparam*r*C*exp(-Dparam*r^2)
        derivative.data[diagind(N+1,N+1)[2:end]] .= -(B*A*sech(B*r)^2)
    end
end
