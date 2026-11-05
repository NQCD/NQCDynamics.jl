using LinearAlgebra: norm

"""
    DiatomicHarmonic(r₀=1.0)

Harmonic interaction between two particles.

Only the first two atoms (the first two columns of `R`) interact. Any additional
atoms are treated as independent and do not contribute to the potential or derivative.
Requires at least 2 atoms.
"""
Parameters.@with_kw struct DiatomicHarmonic <: ClassicalModel
    r₀::Float64 = 1.0
end

NQCModels.ndofs(::DiatomicHarmonic) = 3

# GPU extensions: ext/NQCModelsCUDAExt.jl:54, ext/NQCModelsMetalExt.jl:54
function NQCModels.potential(model::DiatomicHarmonic, R::AbstractMatrix)
    return (norm(R[:,1] .- R[:,2]) - model.r₀)^2 / 2
end

function NQCModels.potential!(model::DiatomicHarmonic, V::Matrix{<:Number}, R::AbstractMatrix)
    V .= NQCModels.potential(model, R)
end

function NQCModels.derivative!(model::DiatomicHarmonic, D::AbstractMatrix, R::AbstractMatrix)
    D .= zero(eltype(D))
    diff = R[:,1] .- R[:,2]
    leng = norm(diff)
    D[:,1] .= (leng - model.r₀) / leng .* diff
    D[:,2] .= .-D[:,1]
end
