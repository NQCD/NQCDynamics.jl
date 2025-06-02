
OrdinaryDiffEq.isfsal(::BCBFull)  = false

struct BCBFullCache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    cayley::Vector{Matrix{uNoUnitsType}}
    halfdt::uNoUnitsType
end

function OrdinaryDiffEq.alg_cache(::BCBFull,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=false)
    halfdt = dt/2
    BCBFullCache(u, uprev, tmp, k, cayley, halfdt)
end

function OrdinaryDiffEq.initialize!(integrator, integrator_cache::BCBFullCache)
    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    NQCCalculators.update_cache!(integrator.p.cache, r)
    DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t)
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, integrator_cache::BCBFullCache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, tmp, cayley, halfdt = integrator_cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)

    rfinal = DynamicsUtils.get_positions(u)
    vfinal = DynamicsUtils.get_velocities(u)
    vtmp = DynamicsUtils.get_velocities(tmp)

    copyto!(rfinal, rprev)

    step_B!(vtmp, vprev, halfdt, k)

    RingPolymerArrays.transform_to_normal_modes!(rfinal, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)
    step_C!(vtmp, rfinal, cayley)
    RingPolymerArrays.transform_from_normal_modes!(rfinal, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)

    NQCCalculators.update_cache!(p.cache, rfinal)
    DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t)

    step_B!(vfinal, vtmp, halfdt, k)
end
