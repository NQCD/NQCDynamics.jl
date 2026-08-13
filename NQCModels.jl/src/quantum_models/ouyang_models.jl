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

function NQCModels.potential(model::OuyangModelOne, R::AbstractMatrix)
    r = R[1]
    Parameters.@unpack A, B, C, D, ΔE, N = model

    V = zeros((N+1, N+1))

    V[1,1] = A*tanh(B*r)
    Δs = range(-ΔE/2, ΔE/2, length=N)
    V[1,2:N+1] .= C*exp(-D*r^2)
    V[diagind(N+1,N+1)[2:end]] .= -(A*tanh(B*r)).+Δs

    return V
end

function NQCModels.potential!(model::OuyangModelOne, V::Hermitian, R::AbstractMatrix)
    r = R[1]
    Parameters.@unpack A, B, C, D, ΔE, N = model

    V[1,1] = A*tanh(B*r)
    Δs = range(-ΔE/2, ΔE/2, length=N)
    V.data[1,2:N+1] .= C*exp(-D*r^2)
    V.data[diagind(N+1,N+1)[2:end]] .= -(A*tanh(B*r)).+Δs
end

function NQCModels.derivative(model::OuyangModelOne, R::AbstractMatrix)
    r = R[1]
    Parameters.@unpack A, B, C, D, N = model

    derivative = Hermitian(zeros((N+1, N+1)))

    derivative.data[1,1] = B*A*sech(B*r)^2
    derivative.data[1,2:N+1] .= -2D*r*C*exp(-D*r^2)
    derivative.data[diagind(N+1,N+1)[2:end]] .= -(B*A*sech(B*r)^2)

    return derivative
end

function NQCModels.derivative!(model::OuyangModelOne, derivative::Hermitian, R::AbstractMatrix)
    r = R[1]
    Parameters.@unpack A, B, C, D, N = model

    derivative.data[1,1] = B*A*sech(B*r)^2
    derivative.data[1,2:N+1] .= -2D*r*C*exp(-D*r^2)
    derivative.data[diagind(N+1,N+1)[2:end]] .= -(B*A*sech(B*r)^2)
end
