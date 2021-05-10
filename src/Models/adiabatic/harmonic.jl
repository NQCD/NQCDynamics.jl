export Harmonic

"""
$(TYPEDEF)

Adiabatic harmonic potential

``V(R) = 1/2 m \\omega^2 R^2``

$(FIELDS)
"""
struct Harmonic <: AdiabaticModel
    "Mass of oscillator"
    mass::Float64
    "Frequency"
    ω::Float64
    "Location of minimum"
    r₀::Float64
end

Harmonic(;mass=1, ω=1, r₀=0) = Harmonic(mass, ω, r₀)

function potential!(model::Harmonic, V::Vector, R::AbstractMatrix)
    V .= sum(0.5 * model.mass * model.ω^2 .* (R .- model.r₀) .^2)
end

function derivative!(model::Harmonic, D::AbstractMatrix, R::AbstractMatrix) 
    D .= model.mass * model.ω^2 .* (R .- model.r₀)
end
