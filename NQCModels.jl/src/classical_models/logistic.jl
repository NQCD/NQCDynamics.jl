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

NQCModels.ndofs(harmonic::Logistic) = 1

function NQCModels.potential(model::Logistic, R::AbstractMatrix)
    r = R[1]
    (;L, k, x₀, c, a) = model
    return L / (1 + exp(-k * (a * r - x₀))) + c
end

function NQCModels.potential!(model::Logistic, V::Matrix{<:Number}, R::AbstractMatrix)
    r = R[1]
    (;L, k, x₀, c, a) = model
    V .= L / (1 + exp(-k * (a * r - x₀))) + c
end

function NQCModels.derivative(model::Logistic, R::AbstractMatrix)
    r = R[1]  
    (;L, k, x₀, a) = model
    return L * a * k * exp(k * (a * r - x₀)) / ((1 + exp(k * ( a * r - x₀)))^2) # analytical derivative of the logistic function
end

function NQCModels.derivative!(model::Logistic, D::AbstractMatrix, R::AbstractMatrix)
    r = R[1]  
    (;L, k, x₀, a) = model
    D .= L * a * k * exp(k * (a * r - x₀)) / ((1 + exp(k * ( a * r - x₀)))^2) # analytical derivative of the logistic function
end