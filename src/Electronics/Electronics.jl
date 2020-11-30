module Electronics

using LinearAlgebra: Hermitian
using ..Models

"""
    ElectronicContainer{T<:AbstractFloat}
    
Container for the model and the model derived electronic quantities.
"""
mutable struct ElectronicContainer{T<:AbstractFloat}
    n_DoF::UInt
    n_states::UInt

    V0::T
    D0::Vector{T}
    diabatic_potential::Hermitian{T, Matrix{T}}
    diabatic_derivative::Vector{Hermitian{T, Matrix{T}}}
    eigenvalues::Vector{T}
    eigenvectors::Matrix{T}
    adiabatic_derivative::Matrix{T}
    nonadiabatic_coupling::Matrix{T}

    function ElectronicContainer{T}(n_states::Integer, n_DoF::Integer) where {T<:AbstractFloat}
        V0 = 0.0
        D0 = zeros(n_DoF)
        V = Hermitian(zeros(n_states, n_states))
        D = [Hermitian(zeros(n_states, n_states)) for _=1:n_DoF]
        eigenvalues = zeros(n_states)
        eigenvectors = zeros(n_states, n_states)
        AD = zeros(n_states, n_states)
        NAC = zeros(n_states, n_states)
        new(n_DoF, n_states, V0, D0, V, D, eigenvalues, eigenvectors, AD, NAC)
    end
end

function ElectronicContainer(n_states::Integer, n_DoF::Integer)
    ElectronicContainer{Float64}(n_states, n_DoF)
end

function calculate_derivative!(model::Model, electronics::ElectronicContainer, R::Vector)
    electronics.D0 .= model.get_D0.(R)
end

end # module