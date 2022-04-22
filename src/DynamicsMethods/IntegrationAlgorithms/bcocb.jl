using UnPack: @unpack
using MuladdMacro: @muladd
using DiffEqBase: @..
using LinearAlgebra: Diagonal
using NQCDynamics: RingPolymers, ndofs
using RingPolymerArrays: RingPolymerArray

StochasticDiffEq.alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::BCOCB) = true

struct BCOCBCache{U,A,K,uEltypeNoUnits,F} <: StochasticDiffEq.StochasticDiffEqMutableCache
    tmp::U
    rtmp::A
    vtmp::A
    k::K
    cayley::Vector{Matrix{uEltypeNoUnits}}
    half::uEltypeNoUnits
    friction::F
end

function StochasticDiffEq.alg_cache(::BCOCB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
    tmp = zero(u)
    rtmp = zero(u.x[1])
    vtmp = zero(rtmp)
    k = zeros(size(rtmp))

    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=true)
    half = uEltypeNoUnits(1//2)

    friction = FrictionCache(p, dt)

    BCOCBCache(tmp, rtmp, vtmp, k, cayley, half, friction)
end

struct MDEFCache{flatV,G,N,C,M,S}
    flatvtmp::flatV
    tmp1::flatV
    tmp2::flatV
    gtmp::G
    noise::N
    c1::C
    c2::C
    sqrtmass::M
    σ::S
end

function FrictionCache(sim::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}, dt) 
    DoFs = ndofs(sim)
    atoms = length(sim.atoms)
    beads = length(sim.beads)

    flatvtmp = zeros(DoFs*atoms)
    tmp1 = zero(flatvtmp)
    tmp2 = zero(flatvtmp)

    gtmp = zeros(DoFs*atoms, DoFs*atoms, beads)
    noise = zero(flatvtmp)

    c1 = Diagonal(zeros(DoFs*atoms, DoFs*atoms))
    c2 = zero(c1)
    
    sqrtmass = 1 ./ sqrt.(repeat(sim.atoms.masses; inner=DoFs))
    σ = zero(sqrtmass)

    MDEFCache(flatvtmp, tmp1, tmp2, gtmp, noise, c1, c2, sqrtmass, σ)
end

struct LangevinCache{C,M,S}
    c1::C
    c2::C
    sqrtmass::M
    σ::S
end

function FrictionCache(sim::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}, dt)
    γ0 = sim.method.γ
    γ = [γ0, 2sqrt.(2sim.beads.normal_mode_springs[2:end])...]
    c1 = exp.(-dt.*γ)
    c2 = sqrt.(1 .- c1.^2)
    sqrtmass = 1 ./ repeat(sqrt.(sim.atoms.masses'), ndofs(sim))
    σ = zero(sqrtmass)

    LangevinCache(c1, c2, sqrtmass, σ)
end

function StochasticDiffEq.initialize!(integrator, cache::BCOCBCache)
    @unpack t,uprev,p = integrator
    v = integrator.uprev.x[1]
    r = integrator.uprev.x[2]

    integrator.f.f1(cache.k,v,r,p,t)
end

@muladd function StochasticDiffEq.perform_step!(integrator,cache::BCOCBCache,f=integrator.f)
    @unpack t,dt,sqdt,uprev,u,p,W = integrator
    @unpack rtmp, vtmp, k, cayley, half, friction = cache

    v1 = uprev.x[1]
    r1 = uprev.x[2]

    step_B!(vtmp, v1, half*dt, k)
    @.. rtmp = r1

    RingPolymerArrays.transform_to_normal_modes!(rtmp, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)

    step_C!(vtmp, rtmp, cayley)

    step_O!(friction, integrator, vtmp, rtmp, t+half*dt)

    step_C!(vtmp, rtmp, cayley)

    RingPolymerArrays.transform_from_normal_modes!(rtmp, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)

    @.. u.x[2] = rtmp

    integrator.f.f1(k,vtmp,u.x[2],p,t+dt)
    step_B!(u.x[1], vtmp, half*dt, k)
end
