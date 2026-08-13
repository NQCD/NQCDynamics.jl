
"""
    Parameters.@with_kw struct Morse{T} <: ClassicalModel

# References

[J. Chem. Phys. 88, 4535 (1988)](https://doi.org/10.1063/1.453761)
"""
Parameters.@with_kw struct Morse{T} <: ClassicalModel
    Dₑ::T = 1.0
    x₀::T = 1.0
    a::T  = 1.0
    m::T = 1.0
end

NQCModels.ndofs(harmonic::Morse) = 1

function NQCModels.potential(model::Morse, R::AbstractMatrix)
    r = R[1]
    (;Dₑ, x₀, a) = model
    return Dₑ * (exp(-a*(r-x₀)) - 1)^2
end

function NQCModels.potential!(model::Morse, V::Matrix{<:Number}, R::AbstractMatrix)
    r = R[1]
    (;Dₑ, x₀, a) = model
    V .= Dₑ * (exp(-a*(r-x₀)) - 1)^2
end

function NQCModels.derivative(model::Morse, R::AbstractMatrix)
    r = R[1] 
    (;Dₑ, x₀, a) = model
    return 2Dₑ * (exp(-a*(r-x₀)) - 1) * -a * exp(-a*(r-x₀))
end

function NQCModels.derivative!(model::Morse, D::AbstractMatrix, R::AbstractMatrix)
    r = R[1] 
    (;Dₑ, x₀, a) = model
    @. D = 2Dₑ * (exp(-a*(r-x₀)) - 1) * -a * exp(-a*(r-x₀))
end

"Eq. 44"
function getω₀(model::Morse)
    (;Dₑ, a, m) = model
    return sqrt(2Dₑ * a^2 / m)
end

"Eq. 36"
function getλ(model::Morse)
    (;Dₑ, a, m) = model
    return sqrt(2m*Dₑ) / a
end

"Eq. 43"
function eigenenergy(model::Morse, n)
    λ = getλ(model)
    n > trunc(Int, λ - 1/2) && throw(DomainError(n, "State $n is unbound!")) # Eq. 40
    ω = getω₀(model)
    harmonic = (n + 1/2)
    Eₙ = harmonic - 1/2λ * harmonic^2
    return Eₙ * ω
end