export DiatomicHarmonic

using LinearAlgebra

"""
    DiatomicHarmonic(r₀=1.0)

Harmonic interaction between two particles.
"""
@with_kw struct DiatomicHarmonic <: AdiabaticModel
    r₀::Float64 = 1.0
end

function potential!(model::DiatomicHarmonic, V::Vector, R::AbstractMatrix)
    V .= (norm(R[:,1] .- R[:,2]) - model.r₀)^2 / 2
end

function derivative!(model::DiatomicHarmonic, D::AbstractMatrix, R::AbstractMatrix) 
    diff = R[:,1] .- R[:,2]
    leng = norm(diff)
    D .= (leng - model.r₀) / leng .* diff
    D[:,2] .*= -1
end
