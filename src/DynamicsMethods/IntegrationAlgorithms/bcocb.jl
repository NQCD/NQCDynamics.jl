using UnPack: @unpack
using MuladdMacro: @muladd
using DiffEqBase: @..
using LinearAlgebra: Diagonal
using NQCDynamics: RingPolymers, ndofs
using RingPolymerArrays: RingPolymerArray
# using OrdinaryDiffEqCore: get_fsalfirstlast

StochasticDiffEq.alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::BCOCB) = true

struct BCOCBCache{U,A,K,uEltypeNoUnits,F} <: StochasticDiffEq.StochasticDiffEqMutableCache
    tmp::U
    rtmp::A
    vtmp::A
    k::K
    cayley::Vector{Matrix{uEltypeNoUnits}}
    halfdt::uEltypeNoUnits
    friction::F
end

OrdinaryDiffEqCore.isfsal(::BCOCB) = false

OrdinaryDiffEqCore.get_fsalfirstlast(cache::BCOCBCache, u::Any) = (nothing, nothing)

function StochasticDiffEq.alg_cache(::BCOCB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
    tmp = zero(u)
    rtmp = zero(u.x[1])
    vtmp = zero(rtmp)
    k = zeros(size(rtmp))

    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=true)
    halfdt = uEltypeNoUnits(dt / 2)

    friction = FrictionCache(p, dt)

    BCOCBCache(tmp, rtmp, vtmp, k, cayley, halfdt, friction)
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
    
    sqrtmass = 1 ./ sqrt.(sim.atoms.masses)
    σ = zero(repeat(sqrtmass; inner=DoFs))

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
    sqrtmass = 1 ./ sqrt.(sim.atoms.masses')
    σ = zero(repeat(sqrtmass, ndofs(sim)))

    LangevinCache(c1, c2, sqrtmass, σ)
end

function StochasticDiffEq.initialize!(integrator, integrator_cache::BCOCBCache)
    (;t, p) = integrator
    v, r = integrator.uprev.x
    integrator.f.f1(integrator_cache.k, v, r, p, t)
end

@muladd function StochasticDiffEq.perform_step!(integrator, integrator_cache::BCOCBCache, f=integrator.f)
    (;t,p) = integrator
    (;vtmp, cayley, halfdt, friction) = integrator_cache

    vprev, rprev = integrator.uprev.x
    acceleration = integrator_cache.k
    v, r, vtmp = OrdinaryDiffEq.OrdinaryDiffEqSymplecticRK.alloc_symp_state(integrator)

    copy!(r, rprev)

    step_B!(vtmp, vprev, halfdt, acceleration)

    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(r, p.beads.transformation)

    step_C!(vtmp, r, cayley)
    step_O!(friction, integrator, vtmp, r, t+halfdt)
    step_C!(vtmp, r, cayley)

    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(r, p.beads.transformation)

    integrator.f.f1(acceleration, vtmp, r, p, t)

    step_B!(v, vtmp, halfdt, acceleration)
end
