module Electronics

using LinearAlgebra: Hermitian
using ..Models

"""
    ElectronicContainer{T<:AbstractFloat}
    
Container for the model and the model derived electronic quantities.
"""
mutable struct ElectronicContainer{T<:AbstractFloat}
    model::Model

    n_DoF::UInt

    V0::T
    D0::Vector{T}
    diabatic_potential::Hermitian{T, Matrix{T}}
    diabatic_derivative::Vector{Hermitian{T, Matrix{T}}}
    eigenvalues::Vector{T}
    eigenvectors::Matrix{T}
    adiabatic_derivative::Matrix{T}
    nonadiabatic_coupling::Matrix{T}

    function ElectronicContainer{T}(model::Model, n_DoF::Integer) where {T<:AbstractFloat}
        V0 = 0.0
        D0 = zeros(n_DoF)
        n_states = model.n_states
        V = Hermitian(zeros(n_states, n_states))
        D = [Hermitian(zeros(n_states, n_states)) for _=1:n_DoF]
        eigenvalues = zeros(n_states)
        eigenvectors = zeros(n_states, n_states)
        AD = zeros(n_states, n_states)
        NAC = zeros(n_states, n_states)
        new(model, n_DoF, V0, D0, V, D, eigenvalues, eigenvectors, AD, NAC)
    end
end

function ElectronicContainer(model::Model, n_DoF::Integer)
    ElectronicContainer{Float64}(model, n_DoF)
end

function calculate_derivative!(electronics::ElectronicContainer, R::Vector)
    electronics.D0 .= electronics.model.get_D0.(R)
end

end # module