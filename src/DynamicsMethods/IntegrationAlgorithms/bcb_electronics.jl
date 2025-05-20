using RingPolymerArrays: RingPolymerArrays

OrdinaryDiffEq.isfsal(::BCBwithTsit5) = false

mutable struct BCBwithTsit5Cache{uType,rType,vType,rateType,uEltypeNoUnits,E} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rType
    vtmp::vType
    k::rateType
    cayley::Vector{Matrix{uEltypeNoUnits}}
    electronic_integrator::E
end

function OrdinaryDiffEq.alg_cache(alg::BCBwithTsit5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
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

function OrdinaryDiffEq.initialize!(integrator, integrator_cache::BCBwithTsit5Cache)

    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    σprev = DynamicsUtils.get_quantum_subsystem(integrator.u)
    NQCCalculators.update_electronics!(integrator.p.cache, r)
    if integrator.p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, σprev)
    elseif integrator.p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, integrator.p.method.state)
    end
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, integrator_cache::BCBwithTsit5Cache, repeat_step=false)
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

    NQCCalculators.update_electronics!(p.cache, rtmp)
    if p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsUtils.acceleration!(k, vtmp, rtmp, p, t, σprev)
    elseif p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsUtils.acceleration!(k, vtmp, rtmp, p, t, p.method.state)
    end
    step_B!(vfinal, vtmp, dt/2, k)

    copyto!(rfinal, rtmp)

    d = NQCCalculators.get_centroid_nonadiabatic_coupling(p.cache, rfinal)
    vals = NQCCalculators.get_centroid_eigen(p.cache, rfinal).values
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
