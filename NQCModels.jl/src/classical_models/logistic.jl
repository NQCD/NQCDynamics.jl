"""
    Parameters.@with_kw struct Logistic{T} <: ClassicalModel

    The sigmoid function is given by:
            L / (1 + exp(-k * (a * x - x0))) + c
        where:
            L is the supremum of the values of the function
            k is the steepness of the function
            x0 is the x value of the sigmoid midpoint
            c is the y value of the sigmoid midpoint
            a is a scaling factor

"""

Parameters.@with_kw struct Logistic{T} <: ClassicalModel
    L::T = 1.0
    k::T = 1.0
    x₀::T = 1.0
    c::T = 1.0
    a::T = 1.0
end

NQCModels.ndofs(::Logistic) = 1

# GPU extensions: ext/NQCModelsCUDAExt.jl:43, ext/NQCModelsMetalExt.jl:43
function NQCModels.potential(model::Logistic, R::AbstractMatrix)
    (;L, k, x₀, c, a) = model
    return sum(L / (1 + exp(-k * (a * x - x₀))) + c for x in R)
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:48, ext/NQCModelsMetalExt.jl:48
function NQCModels.derivative!(model::Logistic, D::AbstractMatrix, R::AbstractMatrix)
    (;L, k, x₀, a) = model
    @. D = L * a * k * exp(k * (a * R - x₀)) / ((1 + exp(k * (a * R - x₀)))^2) # analytical derivative of the logistic function
    return D
end