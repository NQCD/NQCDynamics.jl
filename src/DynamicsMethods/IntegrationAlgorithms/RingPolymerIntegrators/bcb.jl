
using OrdinaryDiffEq
using RingPolymerArrays: RingPolymerArrays
using OrdinaryDiffEqCore: OrdinaryDiffEqCore, get_fsalfirstlast, OrdinaryDiffEqAlgorithm
using StochasticDiffEq: StochasticDiffEq
using OrdinaryDiffEqSymplecticRK


struct BCB <: OrdinaryDiffEqAlgorithm end
struct BCBFull <: OrdinaryDiffEqAlgorithm end

struct BCBwithTsit5{T<:OrdinaryDiffEqAlgorithm} <: OrdinaryDiffEqAlgorithm
    electronic_algorithm::T
end

struct BCBWavefunction <: OrdinaryDiffEqAlgorithm end

mutable struct BCBCache{uType,rateType,uEltypeNoUnits} <:
               OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    fsalfirst::rateType
    halfdt::uEltypeNoUnits
    cayley::Vector{Matrix{uEltypeNoUnits}}
end

OrdinaryDiffEqCore.isfsal(::BCB) = true

OrdinaryDiffEqCore.get_fsalfirstlast(cache::BCBCache, _::Any) = (cache.fsalfirst, cache.k)

function OrdinaryDiffEqCore.alg_cache(::BCB, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        inplace::Val{true}, verbose) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(rate_prototype)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    halfdt = dt / 2
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half = false)
    BCBCache(u, uprev, k, tmp, fsalfirst, halfdt, cayley)
end

function OrdinaryDiffEqSymplecticRK.verify_f2(
    f,
    res,
    p,
    q,
    pa,
    t,
    integrator,
    ::BCBCache,
)
    f(res, p, q, pa, t)
    res == p ? res : OrdinaryDiffEqSymplecticRK.throwex(integrator)
end

function OrdinaryDiffEqCore.initialize!(integrator, integrator_cache::BCBCache)
    integrator.fsalfirst = integrator_cache.fsalfirst
    integrator.fsallast = integrator_cache.k

    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast

    duprev, uprev = integrator.uprev.x
    integrator.f.f1(integrator.k[2].x[1], duprev, uprev, integrator.p, integrator.t)
    OrdinaryDiffEqSymplecticRK.verify_f2(
        integrator.f.f2,
        integrator.k[2].x[2],
        duprev,
        uprev,
        integrator.p,
        integrator.t,
        integrator,
        integrator_cache,
    )
end

@muladd function OrdinaryDiffEqCore.perform_step!(
    integrator,
    integrator_cache::BCBCache,
    repeat_step = false,
)

    @unpack t, dt, p = integrator
    (; cayley, halfdt) = integrator_cache

    vprev, rprev, acceleration =
        OrdinaryDiffEqSymplecticRK.load_symp_state(integrator)
    v, r, vtmp = OrdinaryDiffEqSymplecticRK.alloc_symp_state(integrator)

    copy!(r, rprev)

    step_B!(vtmp, vprev, halfdt, acceleration)

    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(r, p.beads.transformation)

    step_C!(vtmp, r, cayley)

    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(r, p.beads.transformation)

    integrator.f.f1(acceleration, vtmp, r, p, t)

    step_B!(v, vtmp, halfdt, acceleration)

    OrdinaryDiffEqSymplecticRK.store_symp_state!(
        integrator,
        integrator_cache,
        acceleration,
        v,
    )
end

mutable struct BCBwithTsit5Cache{uType,rType,vType,rateType,uEltypeNoUnits,E} <: OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rType
    vtmp::vType
    k::rateType
    cayley::Vector{Matrix{uEltypeNoUnits}}
    electronic_integrator::E
end

OrdinaryDiffEqCore.isfsal(::BCBwithTsit5) = false


OrdinaryDiffEqCore.get_fsalfirstlast(cache::BCBwithTsit5Cache, u::Any) = (nothing, nothing)


function OrdinaryDiffEqCore.alg_cache(alg::BCBwithTsit5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true},verbose) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    electronic_problem = DynamicsMethods.DensityMatrixODEProblem(
        Array(DynamicsUtils.get_quantum_subsystem(u)),
        (0.0, dt),
        NQCModels.nstates(p)
    )
    electronic_integrator = SciMLBase.init(
        electronic_problem, alg.electronic_algorithm;
        save_on=false, save_everystep=false, dt=dt/5, adaptive=false
    )

    tmp = zero(u)
    rtmp = zero(DynamicsUtils.get_positions(u))
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=false)
    BCBwithTsit5Cache(u, uprev, tmp, rtmp, vtmp, k, cayley, electronic_integrator)
end

function OrdinaryDiffEqCore.initialize!(integrator, integrator_cache::BCBwithTsit5Cache)

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

@muladd function OrdinaryDiffEqCore.perform_step!(integrator, integrator_cache::BCBwithTsit5Cache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, rtmp, vtmp, cayley, electronic_integrator = integrator_cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)
    σprev = DynamicsUtils.get_quantum_subsystem(uprev)

    rfinal = DynamicsUtils.get_positions(u)
    vfinal = DynamicsUtils.get_velocities(u)

    copyto!(rtmp, rprev)

    step_B!(vtmp, vprev, dt/2, k)
    RingPolymerArrays.transform_to_normal_modes!(rtmp, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)
    step_C!(vtmp, rtmp, cayley)
    RingPolymerArrays.transform_from_normal_modes!(rtmp, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)

    NQCCalculators.update_cache!(p.cache, rtmp)
    if p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsUtils.acceleration!(k, vtmp, rtmp, p, t, σprev)
    elseif p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsUtils.acceleration!(k, vtmp, rtmp, p, t, p.method.state)
    end
    step_B!(vfinal, vtmp, dt/2, k)

    copyto!(rfinal, rtmp)

    d = NQCCalculators.get_centroid_nonadiabatic_coupling(p.cache, rfinal)
    vals = NQCCalculators.get_centroid_eigen(p.cache, rfinal).w
    DynamicsMethods.update_parameters!(
        electronic_integrator.p,
        vals,
        d,
        RingPolymerArrays.get_centroid(vfinal),
        t+dt
    )

    density_matrix = DynamicsUtils.get_quantum_subsystem(u)
    set_ut!(electronic_integrator, density_matrix, t)
    SciMLBase.step!(electronic_integrator, dt, true)
    copy!(density_matrix, electronic_integrator.u)
end


struct BCBFullCache{uType,rateType,uNoUnitsType} <: OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    k::rateType
    cayley::Vector{Matrix{uNoUnitsType}}
    halfdt::uNoUnitsType
end

OrdinaryDiffEqCore.isfsal(::BCBFull) = false

OrdinaryDiffEqCore.get_fsalfirstlast(cache::BCBFullCache, u::Any) = (nothing, nothing)

function OrdinaryDiffEqCore.alg_cache(::BCBFull,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=false)
    halfdt = dt/2
    BCBFullCache(u, uprev, tmp, k, cayley, halfdt)
end

function OrdinaryDiffEqCore.initialize!(integrator, integrator_cache::BCBFullCache)
    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    NQCCalculators.update_cache!(integrator.p.cache, r)
    DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t)
end

@muladd function OrdinaryDiffEqCore.perform_step!(integrator, integrator_cache::BCBFullCache, repeat_step=false)
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

mutable struct BCBWavefunctionCache{uType,vType,rateType,uEltypeNoUnits} <: OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    vtmp::vType
    k::rateType
    cayley::Vector{Matrix{uEltypeNoUnits}}
end

OrdinaryDiffEqCore.isfsal(::BCBWavefunction) = false


OrdinaryDiffEqCore.get_fsalfirstlast(cache::BCBWavefunctionCache, u::Any) = (nothing, nothing)

function OrdinaryDiffEqCore.alg_cache(::BCBWavefunction,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true},verbose) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=false)
    BCBWavefunctionCache(u, uprev, tmp, vtmp, k, cayley)
end

function OrdinaryDiffEqCore.initialize!(integrator, integrator_cache::BCBWavefunctionCache)
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

@muladd function OrdinaryDiffEqCore.perform_step!(integrator, integrator_cache::BCBWavefunctionCache, repeat_step=false)
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

