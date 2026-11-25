
"""
    struct ErpenbeckThoss{T<:AbstractFloat} <: QuantumModel{Diabatic}

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
struct ErpenbeckThoss{T<:AbstractFloat} <: QuantumModel{Diabatic}
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

    return ErpenbeckThoss(austrip(Γ), morse, D₁, D₂, x₀′, a′, c, V∞, q, ã, x̃, V̄ₖ)
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:352, ext/NQCModelsMetalExt.jl:352
function NQCModels.potential!(model::ErpenbeckThoss, V::Hermitian, R::AbstractMatrix)
    (;morse, c, D₁, D₂, x₀′, a′, V∞, q, ã, x̃, V̄ₖ) = model

    V11 = sum(r -> NQCModels.potential(morse, hcat(r)) + c, R)
    V22 = sum(r -> D₁*exp(-2a′*(r-x₀′)) - D₂*exp(-a′*(r-x₀′)) + V∞, R)
    V12 = sum(r -> V̄ₖ * ((1-q)/2*(1 - tanh((r-x̃)/ã)) + q), R)

    V.data .= [V11 V12; V12 V22]
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:363, ext/NQCModelsMetalExt.jl:363
function NQCModels.derivative!(model::ErpenbeckThoss, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    (;morse, D₁, D₂, x₀′, a′, q, ã, x̃, V̄ₖ) = model
    full_D11 = NQCModels.derivative(morse, R)
    for dof in axes(R, 1)
        for atom in axes(R, 2)
            r = R[dof,atom]
            D11 = full_D11[dof,atom]
            D22 = -2a′*D₁*exp(-2a′*(r-x₀′)) + a′*D₂*exp(-a′*(r-x₀′))
            D12 = -V̄ₖ * (1-q)/2 * sech((r-x̃)/ã)^2 / ã

            D[dof,atom].data[1,1] = D11
            D[dof,atom].data[2,2] = D22
            D[dof,atom].data[1,2] = D12
        end
    end
    return nothing
end

NQCModels.nstates(::ErpenbeckThoss) = 2
NQCModels.ndofs(::ErpenbeckThoss) = 1
