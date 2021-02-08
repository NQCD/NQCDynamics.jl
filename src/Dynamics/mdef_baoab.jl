using UnPack, MuladdMacro
using DiffEqBase: @..
using StochasticDiffEq
using StochasticDiffEq: StochasticDiffEqConstantCache, StochasticDiffEqMutableCache
import StochasticDiffEq: alg_compatible, alg_cache, initialize!, perform_step!

export MDEF_BAOAB

struct MDEF_BAOAB <: StochasticDiffEqAlgorithm end

alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::MDEF_BAOAB) = true

struct MDEF_BAOABConstantCache{uType,uEltypeNoUnits} <: StochasticDiffEqConstantCache
    k::uType
    half::uEltypeNoUnits
end
struct MDEF_BAOABCache{uType,uTypeFlat,uEltypeNoUnits,rateNoiseType,compoundType} <: StochasticDiffEqMutableCache
    utmp::uType
    dutmp::uType
    k::uType
    flatdutmp::uTypeFlat
    gtmp::compoundType
    noise::rateNoiseType
    half::uEltypeNoUnits
    c1::Matrix{uEltypeNoUnits}
    c2::Matrix{uEltypeNoUnits}
end

function alg_cache(alg::MDEF_BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
    k = zero(rate_prototype.x[1])
    MDEF_BAOABConstantCache(k, uEltypeNoUnits(1//2))
end

function alg_cache(alg::MDEF_BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
    dutmp = zero(u.x[1])
    utmp = zero(u.x[2])
    k = zero(rate_prototype.x[1])
    flatdutmp = zero(dutmp[:])

    gtmp = ArrayPartition(zeros(length(utmp), length(utmp)), zeros(length(utmp)))
    noise = zero(rate_prototype.x[1])

    half = uEltypeNoUnits(1//2)
    c1 = zeros(length(utmp), length(utmp))
    c2 = zero(c1)

    MDEF_BAOABCache(utmp, dutmp, k, flatdutmp, gtmp, noise, half, c1, c2)
end

function initialize!(integrator, cache::MDEF_BAOABConstantCache)
    @unpack t,uprev,p = integrator
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    cache.k .= integrator.f.f1(du1,u1,p,t)
end

function initialize!(integrator, cache::MDEF_BAOABCache)
    @unpack t,uprev,p = integrator
    du1 = integrator.uprev.x[1]
    u1 = integrator.uprev.x[2]

    integrator.f.f1(cache.k,du1,u1,p,t)
end

@muladd function perform_step!(integrator,cache::MDEF_BAOABConstantCache,f=integrator.f)
    @unpack t,dt,uprev,p,W = integrator
    @unpack k, half = cache
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    # B
    du2 = du1 + half*dt*k

    # A
    u2 = u1 + half*dt*du2

    # O
    Λ, σ = integrator.g(u2,p,t+dt*half)
    γ, c = eigen(Λ)
    c1 = diagm(exp.(-γ.*dt))
    c2 = diagm(sqrt.(1 .- diag(c1).^2))
    noise = σ.*W.dW
    du3 = c*c1*c'*du2[:] + c*c2*c'*noise[:]
    du3 = reshape(du3, size(du2))

    # A
    u = u2 + half*dt*du3

    # B
    k .= f.f1(du3,u,p,t+dt)
    du = du3 + half*dt*k

    integrator.u = ArrayPartition((du, u))
end

@muladd function perform_step!(integrator,cache::MDEF_BAOABCache,f=integrator.f)
    @unpack t,dt,uprev,u,p,W = integrator
    @unpack utmp, dutmp, k, flatdutmp, gtmp, noise, half, c1, c2 = cache
    du1 = uprev.x[1]
    u1 = uprev.x[2]

    # B
    @.. dutmp = du1 + half*dt*k

    # A
    @.. utmp = u1 + half*dt*dutmp

    # O
    integrator.g(gtmp,utmp,p,t+dt*half)
    Λ = gtmp.x[1]
    σ = gtmp.x[2][1]
    γ, c = eigen(Λ)
    d = diagind(c1)
    for (j,i) in enumerate(diagind(c1))
        c1[i] = exp(-γ[j]*dt)
        c2[i] = sqrt(1 - c1[i]^2)
    end
    @.. noise = σ*W.dW
    flatdutmp .= c*c1*c'*dutmp[:] + c*c2*c'*noise[:]
    dutmp .= reshape(flatdutmp, size(dutmp))

    # A
    @.. u.x[2] = utmp + half*dt*dutmp

    # B
    f.f1(k,dutmp,u.x[2],p,t+dt)
    @.. u.x[1] = dutmp + half*dt*k
end
