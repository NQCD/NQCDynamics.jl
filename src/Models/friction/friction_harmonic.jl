using Random

export FrictionHarmonic

struct FrictionHarmonic <: AdiabaticFrictionModel
    mass::Float64
    ω::Float64
    r₀::Float64
end

FrictionHarmonic(;mass=1, ω=1, r₀=0) = FrictionHarmonic(mass, ω, r₀)

function potential!(model::FrictionHarmonic, V::Vector, R::AbstractMatrix)
    V .= sum(0.5 * model.mass * model.ω^2 .* (R .- model.r₀) .^2)
end

function derivative!(model::FrictionHarmonic, D::AbstractMatrix, R::AbstractMatrix) 
    D .= model.mass * model.ω^2 .* (R .- model.r₀)
end

function friction!(::FrictionHarmonic, F::AbstractMatrix, ::AbstractMatrix)
    randn!(F)
    F .= F'F
    F .= (F + F')/2
end
