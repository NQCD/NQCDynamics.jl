
"""
    struct AndersonHaldane{T<:AbstractFloat} <: QuantumModel

1D two-state diabatic system capable of modelling a Hydrogen atom scattering on a Ge(111) surface.

This model is inspired by ErpenbeckThoss model at J. Chem. Phys. 151, 191101 (2019).

## References
- J. Chem. Phys. 151, 191101 (2019)
"""
struct AndersonHaldane{T<:AbstractFloat} <: QuantumModel
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

    return AndersonHaldane(morse, c, logistic, q, L, x̃₀, Ã, scaledown)
end

function NQCModels.potential(model::AndersonHaldane, r::AbstractMatrix)

    r = r[1] # r is a 1x1 matrix, we only need the first element
    # position-dependent of the coupling function
    (;morse, logistic, scaledown, c) = model
    ϵ₀(x) = (NQCModels.potential(morse, x) + c) * scaledown## U_0 Morse ##
    h(x) = NQCModels.potential(logistic, x)                    ## U_0 Logistic ##
    ϵ₁(x) = ϵ₀(x) + h(x)                                       ## U_1 ##

    (;q, L, x̃₀, Ã) = model
    Vₖ(x) = Ã * ((1-q)/2*(1 - tanh((x-x̃₀)/L)) + q) # coupling energy eq(20) https://journals.aps.org/prb/pdf/10.1103/PhysRevB.97.235452

    V11 = ϵ₀(r)
    V22 = ϵ₁(r)
    V12 = Vₖ(r)
    return Hermitian([V11 V12; V12 V22]) # Diabatic Hamiltonian
end

function NQCModels.potential!(model::AndersonHaldane, V::Hermitian, r::AbstractMatrix)
    r = r[1] # r is a 1x1 matrix, we only need the first element
    # position-dependent of the coupling function
    (;morse, logistic, scaledown, c) = model
    ϵ₀(x) = (NQCModels.potential(morse, x) + c) * scaledown## U_0 Morse ##
    h(x) = NQCModels.potential(logistic, x)                    ## U_0 Logistic ##
    ϵ₁(x) = ϵ₀(x) + h(x)                                       ## U_1 ##

    (;q, L, x̃₀, Ã) = model
    Vₖ(x) = Ã * ((1-q)/2*(1 - tanh((x-x̃₀)/L)) + q) # coupling energy eq(20) https://journals.aps.org/prb/pdf/10.1103/PhysRevB.97.235452

    V11 = ϵ₀(r)
    V22 = ϵ₁(r)
    V12 = Vₖ(r)
    V.data .= [V11 V12; V12 V22]
end

function NQCModels.derivative(model::AndersonHaldane, r::AbstractMatrix)

    r = r[1] # r is a 1x1 matrix, we only need the first element
    # explicit derivative from the .potential above
    (;morse,logistic) = model
    ∂ϵ₀(x) = NQCModels.derivative(morse, x)
    ∂ϵ₁(x) = NQCModels.derivative(morse, x) + NQCModels.derivative(logistic, x)

    (;q, L, x̃₀, Ã) = model
    ∂Vₖ(x) = -Ã * (1-q)/2 * sech((x-x̃₀)/L)^2 / L

    D11 = ∂ϵ₀(r)
    D22 = ∂ϵ₁(r)
    D12 = ∂Vₖ(r)
    return Hermitian([D11 D12; D12 D22])
end

function NQCModels.derivative!(model::AndersonHaldane, D::Hermitian, r::AbstractMatrix)
    r = r[1] # r is a 1x1 matrix, we only need the first element
    # explicit derivative from the .potential above
    (;morse,logistic) = model
    ∂ϵ₀(x) = NQCModels.derivative(morse, x)
    ∂ϵ₁(x) = NQCModels.derivative(morse, x) + NQCModels.derivative(logistic, x)

    (;q, L, x̃₀, Ã) = model
    ∂Vₖ(x) = -Ã * (1-q)/2 * sech((x-x̃₀)/L)^2 / L

    D11 = ∂ϵ₀(r)
    D22 = ∂ϵ₁(r)
    D12 = ∂Vₖ(r)
    D.data .= [D11 D12; D12 D22]
end

NQCModels.nstates(::AndersonHaldane) = 2
NQCModels.ndofs(::AndersonHaldane) = 1
