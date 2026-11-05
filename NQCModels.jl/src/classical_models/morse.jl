
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

NQCModels.ndofs(::Morse) = 1

# GPU extensions: ext/NQCModelsCUDAExt.jl:29, ext/NQCModelsMetalExt.jl:29
function NQCModels.potential(model::Morse{<:Number}, R::AbstractMatrix)
    (;Dₑ, x₀, a) = model
    return sum(Dₑ * (exp(-a*(x-x₀)) - 1)^2 for x in R)
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:35, ext/NQCModelsMetalExt.jl:35
function NQCModels.derivative!(model::Morse{<:Number}, D::AbstractMatrix, R::AbstractMatrix)
    (;Dₑ, x₀, a) = model
    @. D = 2Dₑ * (exp(-a*(R-x₀)) - 1) * -a * exp(-a*(R-x₀))
    return D
end

"Eq. 44"
function getω₀(model::Morse{<:Number})
    (;Dₑ, a, m) = model
    return sqrt(2Dₑ * a^2 / m)
end

"Eq. 36"
function getλ(model::Morse{<:Number})
    (;Dₑ, a, m) = model
    return sqrt(2m*Dₑ) / a
end

"Eq. 43"
function eigenenergy(model::Morse{<:Number}, n)
    λ = getλ(model)
    n > trunc(Int, λ - 1/2) && throw(DomainError(n, "State $n is unbound!")) # Eq. 40
    ω = getω₀(model)
    harmonic = (n + 1/2)
    Eₙ = harmonic - 1/2λ * harmonic^2
    return Eₙ * ω
end
