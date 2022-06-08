using UnPack: @unpack
using MuladdMacro: @muladd
using DiffEqBase: DiffEqBase
using FastBroadcast: @..
using LinearAlgebra: LAPACK, diagm, diag, mul!, diagind
using RecursiveArrayTools: ArrayPartition

using NQCDynamics: get_temperature

StochasticDiffEq.alg_compatible(::DiffEqBase.AbstractSDEProblem,::MDEF_BAOAB) = true

struct MDEF_BAOABConstantCache{uType,uEltypeNoUnits} <: StochasticDiffEq.StochasticDiffEqConstantCache
    k::uType
    half::uEltypeNoUnits
end
struct MDEF_BAOABCache{uType,rType,vType,uTypeFlat,uEltypeNoUnits,rateNoiseType,compoundType} <: StochasticDiffEq.StochasticDiffEqMutableCache
    tmp::uType
    utmp::rType
    dutmp::vType
    k::vType
    flatdutmp::uTypeFlat
    tmp1::uTypeFlat
    tmp2::uTypeFlat
    gtmp::compoundType
    noise::rateNoiseType
    half::uEltypeNoUnits
    c1::Matrix{uEltypeNoUnits}
    c2::Matrix{uEltypeNoUnits}
end

function StochasticDiffEq.alg_cache(::MDEF_BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
    k = zero(rate_prototype.x[1])
    MDEF_BAOABConstantCache(k, uEltypeNoUnits(1//2))
end

function StochasticDiffEq.alg_cache(::MDEF_BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
    tmp = zero(u)
    dutmp = zero(u.x[1])
    utmp = zero(u.x[2])
    k = zero(rate_prototype.x[1])
    flatdutmp = zero(dutmp[:])
    tmp1 = zero(flatdutmp)
    tmp2 = zero(flatdutmp)

    gtmp = zeros(length(utmp), length(utmp))
    noise = zero(rate_prototype.x[1][:])

    half = uEltypeNoUnits(1//2)
    c1 = zeros(length(utmp), length(utmp))
    c2 = zero(c1)

    MDEF_BAOABCache(tmp, utmp, dutmp, k, flatdutmp, tmp1, tmp2, gtmp, noise, half, c1, c2)
end

function StochasticDiffEq.initialize!(integrator, cache::MDEF_BAOABConstantCache)
    @unpack t,uprev,p = integrator
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    cache.k .= integrator.f.f1(du1,u1,p,t)
end

function StochasticDiffEq.initialize!(integrator, cache::MDEF_BAOABCache)
    @unpack t,uprev,p = integrator
    du1 = integrator.uprev.x[1]
    u1 = integrator.uprev.x[2]

    integrator.f.f1(cache.k,du1,u1,p,t)
end

@muladd function StochasticDiffEq.perform_step!(integrator,cache::MDEF_BAOABConstantCache,f=integrator.f)
    @unpack t,dt,sqdt,uprev,p,W = integrator
    @unpack k, half = cache
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    # B
    du2 = du1 + half*dt*k

    # A
    u2 = u1 + half*dt*du2

    # O
    Λ = integrator.g(u2,p,t+dt*half)
    σ = sqrt(get_temperature(p, t+dt*half)) ./ sqrt.(repeat(p.atoms.masses; inner=ndofs(p)))
    γ, c = LAPACK.syev!('V', 'U', Λ) # symmetric eigen
    clamp!(γ, 0, Inf)
    c1 = diagm(exp.(-γ.*dt))
    c2 = diagm(sqrt.(1 .- diag(c1).^2))
    noise = σ.*W.dW[:] / sqdt
    du3 = c*c1*c'*du2[:] + c*c2*c'*noise
    du3 = reshape(du3, size(du2))

    # A
    u = u2 + half*dt*du3

    # B
    k .= f.f1(du3,u,p,t+dt)
    du = du3 + half*dt*k

    integrator.u = ArrayPartition((du, u))
end

@muladd function StochasticDiffEq.perform_step!(integrator,cache::MDEF_BAOABCache,f=integrator.f)
    @unpack t,dt,sqdt,uprev,u,p,W = integrator
    @unpack utmp, dutmp, k, half, gtmp = cache
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    step_B!(dutmp, du1, half*dt, k)

    step_A!(utmp, u1, half*dt, dutmp)

    integrator.g(gtmp,utmp,p,t+dt*half)
    step_O!(cache, integrator)

    step_A!(u.x[2], utmp, half*dt, dutmp)

    integrator.f.f1(k,dutmp,u.x[2],p,t+dt)
    step_B!(u.x[1], dutmp, half*dt, k)
end
