
"""
    struct ErpenbeckThoss{T<:AbstractFloat} <: QuantumModel

1D two-state Quantum system capable of modelling a molecule adsorbed on a metal surface
or a single-molecule junction.

In the two references, all of the parameters are identical except for the particle mass `m`
and the vertical shift `c` applied to the ϵ₀ state.
Both references modify the shift to ensure the quantum ground-state has an energy of 0 eV.
Note that the mass `m` is specified in atomic mass units (amu) **not** atomic units.
If a value for the vertical offset `c` is not explicitly provided whne constructing the model,
it is automatically determined in the constructor from the Morse potential zero-point energy.

# References

- PHYSICAL REVIEW B 97, 235452 (2018)
- J. Chem. Phys. 151, 191101 (2019)
"""
struct ErpenbeckThoss{T<:AbstractFloat} <: QuantumModel
    Γ::T
    morse::ClassicalModels.Morse{T}
    D₁::T
    D₂::T
    x₀′::T
    a′::T
    c::T
    V∞::T
    q::T
    ã::T
    x̃::T
    V̄ₖ::T
end

function ErpenbeckThoss(;
    Γ,

    m   = austrip(1.0u"u"),
    Dₑ  = austrip(3.52u"eV"),
    x₀  = austrip(1.78u"Å"),
    a   = austrip(1.7361u"Å^-1"),
    D₁  = austrip(4.52u"eV"),
    D₂  = austrip(0.79u"eV"),
    x₀′ = x₀,
    a′  = austrip(1.379u"Å^-1"),
    V∞  = austrip(-1.5u"eV"),

    q   = 0.05,
    ã   = austrip(0.5u"Å"),
    x̃   = austrip(3.5u"Å"),
    V̄ₖ  = sqrt(austrip(Γ)/2π),
    c   = nothing
) 
    morse = ClassicalModels.Morse(;Dₑ, x₀, a, m)
    if isnothing(c)
        c = -NQCModels.ClassicalModels.eigenenergy(morse, 0) # Set c to offset zero-point energy
    end

    return ErpenbeckThoss(austrip(Γ), morse, D₁, D₂, x₀′, a′, c, V∞, q, ã, x̃, V̄ₖ)
end

function NQCModels.potential(model::ErpenbeckThoss, R::AbstractMatrix)
    (;morse, c, D₁, D₂, x₀′, a′, V∞) = model

    (;q, ã, x̃, V̄ₖ) = model

    V11 = NQCModels.potential(morse, R) + c
    V22 = D₁*exp(-2a′*(R[1]-x₀′)) - D₂*exp(-a′*(R[1]-x₀′)) + V∞
    V12 = V̄ₖ * ((1-q)/2*(1 - tanh((R[1]-x̃)/ã)) + q)
    return Hermitian([V11 V12; V12 V22])
end

function NQCModels.potential!(model::ErpenbeckThoss, V::Hermitian, R::AbstractMatrix)
    (;morse, c, D₁, D₂, x₀′, a′, V∞) = model

    (;q, ã, x̃, V̄ₖ) = model

    V11 = NQCModels.potential(morse, R) + c
    V22 = D₁*exp(-2a′*(R[1]-x₀′)) - D₂*exp(-a′*(R[1]-x₀′)) + V∞
    V12 = V̄ₖ * ((1-q)/2*(1 - tanh((R[1]-x̃)/ã)) + q)

    V.data[1,1] = V11
    V.data[2,2] = V22
    V.data[1,2] = V12

end

function NQCModels.derivative(model::ErpenbeckThoss, R::AbstractMatrix)
    (;morse, D₁, D₂, x₀′, a′) = model

    (;q, ã, x̃, V̄ₖ) = model

    D11 = NQCModels.derivative(morse, R)
    D22 = -2a′*D₁*exp(-2a′*(R[1]-x₀′)) + a′*D₂*exp(-a′*(R[1]-x₀′))
    D12 = -V̄ₖ * (1-q)/2 * sech((R[1]-x̃)/ã)^2 / ã
    return Hermitian([D11 D12; D12 D22])
end

function NQCModels.derivative!(model::ErpenbeckThoss, D::Matrix{<:Hermitian}, R::AbstractMatrix)
    NQCModels.derivative!(model, D[1,1], R)
end

function NQCModels.derivative!(model::ErpenbeckThoss, D::Hermitian, R::AbstractMatrix)
    (;morse, D₁, D₂, x₀′, a′) = model

    (;q, ã, x̃, V̄ₖ) = model

    D11 = NQCModels.derivative(morse, R)
    D22 = -2a′*D₁*exp(-2a′*(R[1]-x₀′)) + a′*D₂*exp(-a′*(R[1]-x₀′))
    D12 = -V̄ₖ * (1-q)/2 * sech((R[1]-x̃)/ã)^2 / ã
    
    D.data[1,1] = D11
    D.data[2,2] = D22
    D.data[1,2] = D12

    return nothing
end

NQCModels.nstates(::ErpenbeckThoss) = 2
NQCModels.ndofs(::ErpenbeckThoss) = 1
