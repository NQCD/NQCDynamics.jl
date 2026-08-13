using LinearAlgebra: norm

"""
    DiatomicHarmonic(r₀=1.0)

Harmonic interaction between two particles.
"""
Parameters.@with_kw struct DiatomicHarmonic <: ClassicalModel
    r₀::Float64 = 1.0
end

NQCModels.ndofs(::DiatomicHarmonic) = 3

function NQCModels.potential(model::DiatomicHarmonic, R::AbstractMatrix)
    return (norm(R[:,1] .- R[:,2]) - model.r₀)^2 / 2
end

function NQCModels.potential!(model::DiatomicHarmonic, V::Matrix{<:Number}, R::AbstractMatrix)
    V .= (norm(R[:,1] .- R[:,2]) - model.r₀)^2 / 2
end

function NQCModels.derivative!(model::DiatomicHarmonic, D::AbstractMatrix, R::AbstractMatrix) 
    diff = R[:,1] .- R[:,2]
    leng = norm(diff)
    D .= (leng - model.r₀) / leng .* diff
    D[:,2] .*= -1
end
