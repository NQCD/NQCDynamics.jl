using RingPolymerArrays: RingPolymerArrays
using OrdinaryDiffEq.OrdinaryDiffEqCore: get_fsalfirstlast

mutable struct BCBWavefunctionCache{uType,vType,rateType,uEltypeNoUnits} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    vtmp::vType
    k::rateType
    cayley::Vector{Matrix{uEltypeNoUnits}}
end

OrdinaryDiffEq.isfsal(::BCBWavefunction) = false

OrdinaryDiffEq.get_fsalfirstlast(cache::BCBWavefunctionCache, u::Any) = (nothing, nothing)

function OrdinaryDiffEq.alg_cache(::BCBWavefunction,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=false)
    BCBWavefunctionCache(u, uprev, tmp, vtmp, k, cayley)
end

function OrdinaryDiffEq.initialize!(integrator, integrator_cache::BCBWavefunctionCache)
    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    σprev = DynamicsUtils.get_quantum_subsystem(integrator.u)
    NQCCalculators.update_cache!(integrator.p.cache, r)
    if integrator.p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, σprev)
    elseif integrator.p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, integrator.p.method.state)
    end
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, integrator_cache::BCBWavefunctionCache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, vtmp, cayley = integrator_cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)
    σprev = DynamicsUtils.get_quantum_subsystem(uprev)

    rfinal = DynamicsUtils.get_positions(u)
    vfinal = DynamicsUtils.get_velocities(u)
    σfinal = DynamicsUtils.get_quantum_subsystem(u)

    copyto!(rfinal, rprev)

    step_B!(vtmp, vprev, dt/2, k)

    RingPolymerArrays.transform_to_normal_modes!(rfinal, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)
    step_C!(vtmp, rfinal, cayley)
    RingPolymerArrays.transform_from_normal_modes!(rfinal, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)

    NQCCalculators.update_cache!(p.cache, rfinal)
    if p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t, σprev)
    elseif p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t, p.method.state)
    end
    step_B!(vfinal, vtmp, dt/2, k)

    DynamicsUtils.propagate_wavefunction!(σfinal, σprev, vprev, rprev, p, dt)

end
