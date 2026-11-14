module NQCModelsCUDAExt

using NQCModels
using NQCModels: ClassicalModels, QuantumModels, FrictionModels
using CUDA
using LinearAlgebra: Hermitian, diagind
using Parameters

# ============================================================================
# Classical Models - GPU array support
# ============================================================================

# Free model
NQCModels.potential(::ClassicalModels.Free, R::CUDA.CuArray{T}) where {T} = zero(T)

function NQCModels.derivative(model::ClassicalModels.Free, R::CUDA.CuArray)
    return CUDA.zeros(eltype(R), size(R))
end

# Harmonic model
function NQCModels.derivative(model::ClassicalModels.Harmonic, R::CUDA.CuArray)
    return model.m .* model.ω.^2 .* (R .- model.r₀)
end

function NQCModels.potential(model::ClassicalModels.Harmonic, R::CUDA.CuArray)
    return sum(0.5 .* model.m .* model.ω.^2 .* (R .- model.r₀) .^2)
end

# Morse model
function NQCModels.potential(model::ClassicalModels.Morse, R::CUDA.CuArray)
    r = R[1]
    (;Dₑ, x₀, a) = model
    return Dₑ * (exp(-a*(r-x₀)) - 1)^2
end

function NQCModels.derivative(model::ClassicalModels.Morse, R::CUDA.CuArray)
    r = R[1]
    (;Dₑ, x₀, a) = model
    result = CUDA.CuArray([2Dₑ * (exp(-a*(r-x₀)) - 1) * -a * exp(-a*(r-x₀))])
    return result
end

# Logistic model
function NQCModels.potential(model::ClassicalModels.Logistic, R::CUDA.CuArray)
    r = R[1]
    return model.V * (1 + tanh(model.γ * r))
end

function NQCModels.derivative(model::ClassicalModels.Logistic, R::CUDA.CuArray)
    r = R[1]
    return CUDA.CuArray([model.V * model.γ * sech(model.γ * r)^2])
end

# DiatomicHarmonic model
# Only the first two atoms interact; any further atoms remain independent (zero contribution).
function NQCModels.potential(model::ClassicalModels.DiatomicHarmonic, R::CUDA.CuArray)
    diff = R[:, 1] .- R[:, 2]
    r = sqrt(sum(diff .^ 2))
    return (r - model.r₀)^2 / 2
end

function NQCModels.derivative(model::ClassicalModels.DiatomicHarmonic, R::CUDA.CuArray)
    diff = R[:, 1] .- R[:, 2]
    r = sqrt(sum(diff .^ 2))
    D = CUDA.zeros(eltype(R), size(R))
    D[:, 1] .= (r - model.r₀) .* diff ./ r
    D[:, 2] .= .-D[:, 1]
    return D
end

# DarlingHollowayElbow model
function NQCModels.potential(model::ClassicalModels.DarlingHollowayElbow, R::CUDA.CuArray)
    (;a, b, c, d, e, f, g) = model
    x = R[1]
    y = R[2]

    V = a + b*x + c*x^2 + d*y + e*y^2 + f*x*y + g*y^3
    return V
end

function NQCModels.derivative(model::ClassicalModels.DarlingHollowayElbow, R::CUDA.CuArray)
    (;b, c, d, e, f, g) = model
    x = R[1]
    y = R[2]

    dVdx = b + 2c*x + f*y
    dVdy = d + 2e*y + f*x + 3g*y^2

    return CUDA.CuArray([dVdx, dVdy])
end

# ============================================================================
# Quantum Models - GPU array support
# ============================================================================

# Helper function for creating GPU-backed matrices
function gpu_hermitian_matrix(model::QuantumModels.QuantumModel, ::Type{T}) where T<:CUDA.CuArray
    n = NQCModels.nstates(model)
    return CUDA.zeros(eltype(T), n, n)
end

function NQCModels.zero_derivative(model::QuantumModels.QuantumModel, R::CUDA.CuArray)
    # For quantum models, return array of regular matrices (not Hermitian wrapper for GPU arrays)
    return [CUDA.zeros(eltype(R), NQCModels.nstates(model), NQCModels.nstates(model)) for _ in CartesianIndices(R)]
end

# TullyModelOne
function NQCModels.potential(model::QuantumModels.TullyModelOne, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d = model
    q = R[1]
    if q > 0
        V11 = a * (1 - exp(-b*q))
    else
        V11 = -a * (1 - exp(b*q))
    end
    V22 = -V11
    V12 = c * exp(-d*q^2)
    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.TullyModelOne, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d = model
    q = R[1]
    D11 = a * b * exp(-b * abs(q))
    D22 = -D11
    D12 = -2 * c * d * q * exp(-d*q^2)
    return CUDA.CuArray([D11 D12; D12 D22])
end

# TullyModelTwo
function NQCModels.potential(model::QuantumModels.TullyModelTwo, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d, e = model
    q = R[1]
    V11 = 0
    V22 = -a*exp(-b*q^2) + e
    V12 = c * exp(-d*q^2)
    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.TullyModelTwo, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d = model
    q = R[1]
    D11 = 0
    D22 = 2*a*b*q*exp(-b*q^2)
    D12 = -2*c*d*q*exp(-d*q^2)
    return CUDA.CuArray([D11 D12; D12 D22])
end

# TullyModelThree
function NQCModels.potential(model::QuantumModels.TullyModelThree, R::CUDA.CuArray)
    Parameters.@unpack a, b, c = model
    q = R[1]
    V11 = a
    V22 = -a
    if q > 0
        V12 = b * (2 - exp(-c*q))
    else
        V12 = b * exp(c*q)
    end
    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.TullyModelThree, R::CUDA.CuArray)
    Parameters.@unpack a, b, c = model
    q = R[1]
    D11 = 0
    D22 = 0
    if q > 0
        D12 = b * c * exp(-c*q)
    else
        D12 = b * c * exp(c*q)
    end
    return CUDA.CuArray([D11 D12; D12 D22])
end

# DoubleWell
function NQCModels.potential(model::QuantumModels.DoubleWell, R::CUDA.CuArray)
    r = R[1]

    V₀ = 0.5 * model.mass * model.ω^2 * r^2
    v = sqrt(2)*model.γ*r
    V11 = V₀ + v
    V22 = V₀ - v
    V12 = model.Δ/2

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.DoubleWell, R::CUDA.CuArray)
    r = R[1]

    D₀ = model.mass * model.ω^2 * r
    v = sqrt(2)*model.γ
    D11 = D₀ + v
    D22 = D₀ - v

    return CUDA.CuArray([D11 0; 0 D22])
end

# AnanthModelOne
function NQCModels.potential(model::QuantumModels.AnanthModelOne, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d = model
    q = R[1]

    V11 = a * tanh(b*q)
    V22 = -a * tanh(b*q)
    V12 = c * exp(-d*q^2)

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.AnanthModelOne, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d = model
    q = R[1]

    D11 = a * b * (1 - tanh(b*q)^2)
    D22 = -D11
    D12 = -2 * c * d * q * exp(-d*q^2)

    return CUDA.CuArray([D11 D12; D12 D22])
end

# AnanthModelTwo
function NQCModels.potential(model::QuantumModels.AnanthModelTwo, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d, e, f = model
    q = R[1]

    V11 = a * tanh(f*q)
    V22 = -b * tanh(f*q)
    V12 = c * exp(-d*(q+e)^2)

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.AnanthModelTwo, R::CUDA.CuArray)
    Parameters.@unpack a, b, c, d, e, f = model
    q = R[1]

    D11 = a * f * (1 - tanh(f*q)^2)
    D22 = -b * f * (1 - tanh(f*q)^2)
    D12 = -2 * c * d * (q+e) * exp(-d*(q+e)^2)

    return CUDA.CuArray([D11 D12; D12 D22])
end

# MiaoSubotnik
function NQCModels.potential(model::QuantumModels.MiaoSubotnik, R::CUDA.CuArray)
    (;m, ω, g, ΔG, coupling) = model
    x = R[1]

    V11 = (m*ω^2*x^2)/2
    V22 = (m*ω^2*(x-g)^2)/2 + ΔG
    V12 = coupling

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.MiaoSubotnik, R::CUDA.CuArray)
    (;m, ω, g) = model
    x = R[1]

    D11 = m*ω^2*x
    D22 = m*ω^2*(x-g)

    return CUDA.CuArray([D11 0; 0 D22])
end

# GatesHollowayElbow
function NQCModels.potential(model::QuantumModels.GatesHollowayElbow, R::CUDA.CuArray)
    Parameters.@unpack λ₁, λ₂, z₀, x₀, α, d, z12, c, γ = model

    x = R[1]
    z = R[2]

    V11 = d*(1-exp(-α*x))^2 + exp(-λ₁*(z+z₀))
    V22 = d*(1-exp(-α*z))^2 + exp(-λ₂*(x+x₀))
    V12 = c * exp(-γ*(z-z12))

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.GatesHollowayElbow, R::CUDA.CuArray)
    Parameters.@unpack λ₁, λ₂, z₀, x₀, α, d, z12, c, γ = model

    x = R[1, 1]
    z = R[2, 1]

    # x derivative
    D11_x = 2α*d*(1-exp(-α*x))*exp(-α*x)
    D22_x = -λ₂*exp(-λ₂*(x+x₀))
    D12_x = 0

    # z derivative
    D11_z = -λ₁*exp(-λ₁*(z+z₀))
    D22_z = 2α*d*(1-exp(-α*z))*exp(-α*z)
    D12_z = c * -γ*exp(-γ*(z-z12))

    D_x = CUDA.CuArray([D11_x D12_x; D12_x D22_x])
    D_z = CUDA.CuArray([D11_z D12_z; D12_z D22_z])

    return [D_x, D_z]
end

# ThreeStateMorse
function NQCModels.potential(model::QuantumModels.ThreeStateMorse, R::CUDA.CuArray)
    TSM_V_ii(x, d, α, r, c) = d * (1 - exp(-α*(x-r)))^2 + c
    TSM_V_ij(x, a, α, r) = a * exp(-α*(x-r)^2)

    r = R[1]

    V11 = TSM_V_ii(r, model.d1, model.α1, model.r1, model.c1)
    V22 = TSM_V_ii(r, model.d2, model.α2, model.r2, model.c2)
    V33 = TSM_V_ii(r, model.d3, model.α3, model.r3, model.c3)

    V12 = TSM_V_ij(r, model.a12, model.α12, model.r12)
    V13 = TSM_V_ij(r, model.a13, model.α13, model.r13)
    V23 = TSM_V_ij(r, model.a23, model.α23, model.r23)

    return CUDA.CuArray([V11 V12 V13; V12 V22 V23; V13 V23 V33])
end

function NQCModels.derivative(model::QuantumModels.ThreeStateMorse, R::CUDA.CuArray)
    function TSM_D_ii(x, d, α, r)
        ex = exp(-α*(x-r))
        return 2 * d * α * (ex - ex^2)
    end
    TSM_D_ij(x, a, α, r) = -2 * a * α * (x-r) * exp(-α*(x-r)^2)

    r = R[1]

    D11 = TSM_D_ii(r, model.d1, model.α1, model.r1)
    D22 = TSM_D_ii(r, model.d2, model.α2, model.r2)
    D33 = TSM_D_ii(r, model.d3, model.α3, model.r3)

    D12 = TSM_D_ij(r, model.a12, model.α12, model.r12)
    D13 = TSM_D_ij(r, model.a13, model.α13, model.r13)
    D23 = TSM_D_ij(r, model.a23, model.α23, model.r23)

    return CUDA.CuArray([D11 D12 D13; D12 D22 D23; D13 D23 D33])
end

# ErpenbeckThoss
function NQCModels.potential(model::QuantumModels.ErpenbeckThoss, R::CUDA.CuArray)
    (;morse, c, D₁, D₂, x₀′, a′, V∞) = model
    (;q, ã, x̃, V̄ₖ) = model

    V11 = NQCModels.potential(morse, R) + c
    V22 = D₁*exp(-2a′*(R[1]-x₀′)) - D₂*exp(-a′*(R[1]-x₀′)) + V∞
    V12 = V̄ₖ * ((1-q)/2*(1 - tanh((R[1]-x̃)/ã)) + q)

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.ErpenbeckThoss, R::CUDA.CuArray)
    (;morse, D₁, D₂, x₀′, a′) = model
    (;q, ã, x̃, V̄ₖ) = model

    D11 = NQCModels.derivative(morse, R)
    D22 = -2a′*D₁*exp(-2a′*(R[1]-x₀′)) + a′*D₂*exp(-a′*(R[1]-x₀′))
    D12 = -V̄ₖ * (1-q)/2 * sech((R[1]-x̃)/ã)^2 / ã

    return CUDA.CuArray([D11 D12; D12 D22])
end

# AndersonHaldane
function NQCModels.potential(model::QuantumModels.AndersonHaldane, R::CUDA.CuArray)
    r = R[1]
    (;morse, logistic, scaledown, c) = model
    ϵ₀ = (NQCModels.potential(morse, R) + c) * scaledown
    h = NQCModels.potential(logistic, R)
    ϵ₁ = ϵ₀ + h

    (;q, L, x̃₀, Ã) = model
    Vₖ = Ã * ((1-q)/2*(1 - tanh((r-x̃₀)/L)) + q)

    V11 = ϵ₀
    V22 = ϵ₁
    V12 = Vₖ

    return CUDA.CuArray([V11 V12; V12 V22])
end

function NQCModels.derivative(model::QuantumModels.AndersonHaldane, R::CUDA.CuArray)
    r = R[1]
    (;morse, logistic) = model
    ∂ϵ₀ = NQCModels.derivative(morse, R)
    ∂ϵ₁ = NQCModels.derivative(morse, R) + NQCModels.derivative(logistic, R)

    (;q, L, x̃₀, Ã) = model
    ∂Vₖ = -Ã * (1-q)/2 * sech((r-x̃₀)/L)^2 / L

    D11 = ∂ϵ₀
    D22 = ∂ϵ₁
    D12 = ∂Vₖ

    return CUDA.CuArray([D11 D12; D12 D22])
end

# OuyangModelOne
function NQCModels.potential(model::QuantumModels.OuyangModelOne, R::CUDA.CuArray)
    r = R[1]
    Parameters.@unpack A, B, C, D, ΔE, N = model

    V = CUDA.zeros((N+1, N+1))

    V[1,1] = A*tanh(B*r)
    Δs = range(-ΔE/2, ΔE/2, length=N)
    V[1,2:N+1] .= C*exp(-D*r^2)
    V[diagind(N+1,N+1)[2:end]] .= -(A*tanh(B*r)).+Δs

    return V
end

function NQCModels.derivative(model::QuantumModels.OuyangModelOne, R::CUDA.CuArray)
    r = R[1]
    Parameters.@unpack A, B, C, D, N = model

    derivative = CUDA.zeros((N+1, N+1))

    derivative[1,1] = B*A*sech(B*r)^2
    derivative[1,2:N+1] .= -2D*r*C*exp(-D*r^2)
    derivative[diagind(N+1,N+1)[2:end]] .= -(B*A*sech(B*r)^2)

    return derivative
end

# Scattering1D
function NQCModels.potential(model::QuantumModels.Scattering1D, R::CUDA.CuArray)
    r = R[1]

    V₀ = model.D*(exp(-2model.a*r) - 2*exp(-model.a*r))
    γᵣ = model.B*exp(-model.a*r^2)

    V = CUDA.zeros((model.n_states, model.n_states))

    V[diagind(V,0)] .= model.α + V₀ #metal state energies
    V[diagind(V,1)] .= model.β #metal-metal couplings between neighbouring states

    V[end-1, end] = 0.0 #periodic coupling between first and last metal state
    V[2,end] = model.β

    V[1,1] = V₀ # Occupied molecule state
    V[1,2] = γᵣ # molecule metal coupling

    return V
end

function NQCModels.derivative(model::QuantumModels.Scattering1D, R::CUDA.CuArray)
    r = R[1]

    D₀ = 2*model.D*model.a*(exp(-model.a*r)-exp(-2*model.a*r))
    dγᵣ = -2model.a*r*model.B*exp(-model.a*r^2)

    D = CUDA.zeros((model.n_states, model.n_states))

    D[1,2] = dγᵣ
    D[diagind(D,0)] .= D₀

    return D
end

# ============================================================================
# Friction Models - GPU array support
# ============================================================================

function FrictionModels.zero_friction(::FrictionModels.Model, R::CUDA.CuArray)
    return CUDA.zeros(eltype(R), length(R), length(R))
end

function FrictionModels.friction(model::FrictionModels.ConstantFriction, R::CUDA.CuArray)
    F = CUDA.zeros(eltype(R), size(R, 1) * size(R, 2), size(R, 1) * size(R, 2))
    if model.γ isa Number
        F .= model.γ
    elseif model.γ isa CUDA.CuArray
        F .= model.γ
    else
        F .= CUDA.CuArray(model.γ)
    end
    return F
end

end # module
