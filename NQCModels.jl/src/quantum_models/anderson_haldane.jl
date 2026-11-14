
"""
    struct AndersonHaldane{T<:AbstractFloat} <: QuantumModel{Diabatic}

1D two-state diabatic system capable of modelling a Hydrogen atom scattering on a Ge(111) surface.

This model is inspired by ErpenbeckThoss model at J. Chem. Phys. 151, 191101 (2019).

## References
- J. Chem. Phys. 151, 191101 (2019)
"""
struct AndersonHaldane{T<:AbstractFloat} <: QuantumModel{Diabatic}
    #Γ::T
    # Morse Potential
    morse::ClassicalModels.Morse{T}
    c::T

    # Logistic Repulsive Curve
    logistic::ClassicalModels.Logistic{T}

    # Coupling tanh Function
    q::T
    L::T
    x̃₀::T
    Ã::T
    scaledown::T
end

function AndersonHaldane(;
    #Γ = 2.18305,

    # Morse Potential
    m   = austrip(1.0u"u"),
    Dₑ  = austrip(0.0502596u"eV"),
    x₀  = austrip(2.07276u"Å"),
    a   = austrip(2.39727u"Å^-1"),
    c  = austrip(2.67532u"eV"),

    # Logistic for h(x)
    D₁ = austrip(4.351u"eV"),
    a′ = austrip(3.9796u"Å^-1"),
    x₀′ = austrip(2.2798u"Å"),
    c′ = austrip(-0.33513u"eV"),
    b = 1.02971,

    # Coupling Function
    Ã  = austrip(2.28499u"eV"),
    L   = austrip(1.10122u"Å"),
    x̃₀   = austrip(2.0589u"Å"),
    q   = 2.26384e-16,

    scaledown = 1.0,

)
    morse = ClassicalModels.Morse(;Dₑ, x₀, a, m)
    logistic = ClassicalModels.Logistic(D₁, a′, x₀′, c′, b)
    #c = -NQCModels.AdiabaticModels.eigenenergy(morse, 0) # Set c to offset zero-point energy

    return AndersonHaldane(morse, c, logistic, q, L, x̃₀, Ã, scaledown)
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:375, ext/NQCModelsMetalExt.jl:375
function NQCModels.potential!(model::AndersonHaldane, V::Hermitian, R::AbstractMatrix)
    (;morse, logistic, scaledown, c, q, L, x̃₀, Ã) = model

    ϵ₀(r) = (NQCModels.potential(morse, hcat(r)) + c) * scaledown        ## U_0 Morse ##
    V11 = sum(ϵ₀, R)
    V22 = sum(r -> ϵ₀(r) + NQCModels.potential(logistic, hcat(r)), R)    ## U_1 = U_0 + U_0 Logistic ##
    # coupling energy eq(20) https://journals.aps.org/prb/pdf/10.1103/PhysRevB.97.235452
    V12 = sum(r -> Ã * ((1-q)/2*(1 - tanh((r-x̃₀)/L)) + q), R)

    V.data .= [V11 V12; V12 V22] # Diabatic Hamiltonian
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:392, ext/NQCModelsMetalExt.jl:392
function NQCModels.derivative!(model::AndersonHaldane, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    (;morse, logistic, q, L, x̃₀, Ã) = model

    ∂ϵ₀_full = NQCModels.derivative(morse, R)
    ∂ϵ₁_full = ∂ϵ₀_full + NQCModels.derivative(logistic, R)
    for dof in axes(R, 1)
        for atom in axes(R, 2)
            r = R[dof,atom]
            ∂ϵ₀ = ∂ϵ₀_full[dof,atom]
            ∂ϵ₁ = ∂ϵ₁_full[dof,atom]
            ∂Vₖ = -Ã * (1-q)/2 * sech((r-x̃₀)/L)^2 / L

            D11 = ∂ϵ₀
            D22 = ∂ϵ₁
            D12 = ∂Vₖ
            D[dof,atom].data .= [D11 D12; D12 D22]
        end
    end
end

NQCModels.nstates(::AndersonHaldane) = 2
NQCModels.ndofs(::AndersonHaldane) = 1
