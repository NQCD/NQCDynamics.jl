# The aim of this file is to build a generic universal integrator for all our
# problems by following the BAOAB framework.
using OrdinaryDiffEqCore: OrdinaryDiffEqCore, OrdinaryDiffEqAlgorithm,
                          OrdinaryDiffEqMutableCache, get_fsalfirstlast,
                          update_coefficients!, alg_cache, @cache
using StochasticDiffEq
using StochasticDiffEq: StochasticDiffEqCore
using SciMLBase: SciMLBase, set_ut!
using DiffEqBase
using MuladdMacro: @muladd
using UnPack: @unpack
using FastBroadcast: @..
using NQCDynamics: DynamicsUtils, get_temperature
using .DynamicsUtils: acceleration!, get_positions, get_velocities, get_quantum_subsystem
using NQCCalculators
using LinearAlgebra: LAPACK, diagm, diag, mul!, diagind
using RecursiveArrayTools: ArrayPartition

# ---------------------------------------------------------------------------
# Noise-process tags: dispatch which concrete algorithm/cache pairing is built
# ---------------------------------------------------------------------------
abstract type AbstractXProcess end
struct NoiseCoupled <: AbstractXProcess end
struct NoiseFree <: AbstractXProcess end

# ---------------------------------------------------------------------------
# Algorithm types
#
# These MUST be separate concrete types (not one type parameterized over T)
# because the problem/solver pairing check in DiffEqBase dispatches on the
# *abstract supertype* of the algorithm. A single BAXAB{T} <: StochasticDiffEqAlgorithm
# would always be rejected for ODEProblems regardless of T.
# ---------------------------------------------------------------------------
struct BAXABSDE{T} <: StochasticDiffEq.StochasticDiffEqAlgorithm
    X::Function
    friction::T

end

abstract type FrictionType end
struct ConstantFriction <: FrictionType
    γ::AbstractFloat
end

struct VariableFriction <: FrictionType 
    γ::Nothing
end 

struct BAXABODE <: OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm
    X::Function
end

"""
    BAXAB(process, X, γ=nothing)

Dispatching constructor: returns a `BAXABSDE` for `NoiseCoupled()` (paired with
SDEProblem) or a `BAXABODE` for `NoiseFree()` (paired with ODEProblem). Callers
use `BAXAB(...)` exactly as before; only the returned concrete type changes.
"""
BAXAB(::NoiseCoupled, X::Function, γ::FrictionType=VariableFriction(nothing)) = BAXABSDE(X, γ)
BAXAB(::NoiseFree,    X::Function) = BAXABODE(X)

const AnyBAXAB = Union{BAXABSDE,BAXABODE}

OrdinaryDiffEqCore.isfsal(::AnyBAXAB) = false
alg_order(::AnyBAXAB) = 2

StochasticDiffEqCore.alg_compatible(::DiffEqBase.AbstractSDEProblem, ::BAXABSDE) = true
StochasticDiffEqCore.alg_compatible(::DiffEqBase.AbstractSDEProblem, ::BAXABODE) = false

# ---------------------------------------------------------------------------
# Caches
#
# Each cache only holds what its corresponding step_X! function actually uses.
# BAXABSDECache: full noise/friction machinery for noisedriven_step_X!.
# BAXABODECache: minimal fields for electrondriven_step_X!.
# ---------------------------------------------------------------------------
@cache mutable struct BAXABSDECache{uType,uTypeFlat,uEltypeNoUnits,rateNoiseType,compoundType,rateType} <: OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    tmp::uType
    k::rateType
    flatdutmp::uTypeFlat
    tmp1::uTypeFlat
    tmp2::uTypeFlat
    gtmp::compoundType
    noise::rateNoiseType
    half::uEltypeNoUnits
    c1::Matrix{uEltypeNoUnits}
    c2::Matrix{uEltypeNoUnits}
end

@cache mutable struct BAXABODECache{uType,uEltypeNoUnits,rateType} <: OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    tmp::uType
    k::rateType
    half::uEltypeNoUnits
end

function StochasticDiffEqCore.alg_cache(alg::BAXABSDE, prob, u, ΔW, ΔZ, p,
        rate_prototype, noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, f, t, dt, ::Type{Val{true}}, verbose) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    k = zero(DynamicsUtils.get_velocities(rate_prototype))

    utmp = DynamicsUtils.get_positions(tmp)
    flatdutmp = zero(vec(DynamicsUtils.get_velocities(tmp)))
    tmp1 = zero(flatdutmp)
    tmp2 = zero(flatdutmp)

    n = length(utmp)
    
    noise = zero(vec(DynamicsUtils.get_velocities(rate_prototype)))

    half = uEltypeNoUnits(1//2)

    gtmp, c1, c2 = alg.friction(f.g, alg.friction.γ, n, dt)

    BAXABSDECache(tmp, k, flatdutmp, tmp1, tmp2, gtmp, noise, half, c1, c2)
end

function (::VariableFriction)(g, γ,n, dt)
    gtmp = zeros(n, n)
    c1 = zeros(n, n)
    c2 = zeros(n, n)
    return gtmp, c1, c2
end

function (::ConstantFriction)(g, γ, n, dt)
    gtmp = zeros(n, n)
    
    c1 = zeros(n, n)
    c2 = zeros(n, n)
    for (j,i) in enumerate(diagind(c1))
        c1[i] = exp(-γ*dt)
        c2[i] = sqrt(1 - c1[i]^2)
    end 
    return gtmp, c1, c2
end


function OrdinaryDiffEqCore.alg_cache(::BAXABODE, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        inplace::Val{true}, verbose) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    k = zero(DynamicsUtils.get_velocities(rate_prototype))
    half = uEltypeNoUnits(1//2)
    BAXABODECache(tmp, k, half)
end

OrdinaryDiffEqCore.get_fsalfirstlast(cache::BAXABSDECache, rate_prototype) = (nothing, nothing)
OrdinaryDiffEqCore.get_fsalfirstlast(cache::BAXABODECache, rate_prototype) = (nothing, nothing)

# ---------------------------------------------------------------------------
# initialize! — package-specific entry points sharing one implementation
# ---------------------------------------------------------------------------
function StochasticDiffEq.initialize!(integrator, cache::BAXABSDECache)
    _baxab_initialize!(integrator, cache)
end

function OrdinaryDiffEqCore.initialize!(integrator, cache::BAXABODECache)
    _baxab_initialize!(integrator, cache)
end

function _baxab_initialize!(integrator, cache)
    @unpack t, uprev, p = integrator
    NQCCalculators.update_cache!(p.cache, DynamicsUtils.get_positions(uprev))
    integrator.f(cache.k, uprev, p, t)
end

# ---------------------------------------------------------------------------
# perform_step! — the BAOAB skeleton only ever touches tmp/k/half, which both
# caches provide. The alg-specific X-step (noise vs electron driven) is where
# the cache-specific fields get used.
# ---------------------------------------------------------------------------
function StochasticDiffEq.perform_step!(integrator, cache::BAXABSDECache, f=integrator.f)
    _baxab_perform_step!(integrator, cache, f)
end

function OrdinaryDiffEqCore.perform_step!(integrator, cache::BAXABODECache, f=integrator.f)
    _baxab_perform_step!(integrator, cache, f)
end

function _baxab_perform_step!(integrator, cache, f)
    @unpack t, dt, uprev, u, f, p = integrator
    @unpack tmp, k, half = cache

    du1 = DynamicsUtils.get_velocities(uprev)
    u1 = DynamicsUtils.get_positions(uprev)
    dutmp = DynamicsUtils.get_velocities(tmp)
    utmp = DynamicsUtils.get_positions(tmp)

    step_AB!(dutmp, du1, half*dt, k)      # B-step
    step_AB!(utmp, u1, half*dt, dutmp)    # A-step
    NQCCalculators.update_cache!(p.cache, utmp)

    integrator.alg.X(integrator, cache)   # X-step (noise- or electron-driven)

    step_AB!(DynamicsUtils.get_positions(u), utmp, half*dt, dutmp)  # A-step

    copyto!(DynamicsUtils.get_velocities(u), dutmp)
    NQCCalculators.update_cache!(p.cache, utmp)
    f.f(k, u, p, t + dt)          # ODE/SDE acceleration function
    
    step_AB!(DynamicsUtils.get_velocities(u), dutmp, half*dt, k)    # B-step
end

# ---------------------------------------------------------------------------
# X-steps: each uses only fields present on its matching cache type.
# ---------------------------------------------------------------------------
function noisedriven_step_X!(integrator, cache::BAXABSDECache)
    @unpack t, dt, uprev, p, W, f = integrator
    @unpack k, tmp, half, gtmp, flatdutmp, tmp1, tmp2, noise, c1, c2 = cache

    r = DynamicsUtils.get_positions(tmp)
    v = DynamicsUtils.get_velocities(tmp)

    #f.f(k, tmp, p, t)                 # acceleration function
    f.g(gtmp, r, p, t + dt*half)       # friction function
    Λ = gtmp
    σ = repeat(sqrt.(get_temperature(p, t + dt*half) ./ p.atoms.masses); inner=ndofs(p))

    @.. noise = σ*W.dW[:] / sqrt(dt)

    γ, c = LAPACK.syev!('V', 'U', Λ)   # symmetric eigendecomposition
    clamp!(γ, 0, Inf)
    for (j, i) in enumerate(diagind(c1))
        c1[i] = exp(-γ[j]*dt)
        c2[i] = sqrt(1 - c1[i]^2)
    end

    copyto!(flatdutmp, v)
    mul!(tmp1, transpose(c), flatdutmp)
    mul!(tmp2, c1, tmp1)
    mul!(flatdutmp, c, tmp2)

    mul!(tmp1, transpose(c), noise)
    mul!(tmp2, c2, tmp1)
    mul!(tmp1, c, tmp2)

    @.. flatdutmp += tmp1
    copyto!(DynamicsUtils.get_velocities(tmp), flatdutmp)
end

function electrondriven_step_X!(integrator, cache::BAXABODECache)
    @unpack t, dt, uprev, u, p, f = integrator
    @unpack k, tmp = cache

    r = DynamicsUtils.get_positions(tmp)
    v = DynamicsUtils.get_velocities(tmp)

    σprev = DynamicsUtils.get_quantum_subsystem(uprev)
    σfinal = DynamicsUtils.get_quantum_subsystem(u)

    DynamicsUtils.propagate_electrons!(σfinal, σprev, v, r, p, dt)
end

function step_AB!(x2, x1, dt, dx)
    @.. broadcast=false x2 = muladd(dt, dx, x1)
end