using LinearAlgebra: diagind

export FreeConstantFriction

struct FreeConstantFriction{T} <: AdiabaticFrictionModel
    γ::T
end

potential!(::FreeConstantFriction, out::Vector, ::AbstractMatrix) = out .= 0
derivative!(::FreeConstantFriction, out::AbstractMatrix, ::AbstractMatrix) = out .= 0

function friction!(model::FreeConstantFriction, F::AbstractMatrix, ::AbstractMatrix)
    F[diagind(F)] .= model.γ
end
