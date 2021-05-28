export Harmonic

"""
    Harmonic(m=1.0, ω=1.0, r₀=0.0)

Adiabatic harmonic potential. ``V(x) = mω^2(x-r₀)^2 / 2``

```jldoctest
julia> using Symbolics;

julia> @variables x, m, ω, r₀;

julia> model = Models.Harmonic(m=m, ω=ω, r₀=r₀);

julia> Models.potential(model, hcat(x))
1-element Vector{Num}:
 0.5m*(ω^2)*((x - r₀)^2)

julia> Models.derivative(model, hcat(x))
1×1 Matrix{Num}:
 m*(x - r₀)*(ω^2)
```
"""
@with_kw struct Harmonic{M,W,R} <: AdiabaticModel
    m::M = 1.0
    ω::W = 1.0
    r₀::R = 0.0
end

function potential!(model::Harmonic, V::Vector, R::AbstractMatrix)
    V .= sum(0.5 * model.m* model.ω^2 .* (R .- model.r₀) .^2)
end

function derivative!(model::Harmonic, D::AbstractMatrix, R::AbstractMatrix) 
    D .= model.m* model.ω^2 .* (R .- model.r₀)
end
