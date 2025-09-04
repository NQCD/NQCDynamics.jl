# using OrdinaryDiffEq.OrdinaryDiffEqCore: get_fsalfirstlast

mutable struct BABwithTsit5Cache{uType,vType,rateType,E} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    vtmp::vType
    k::rateType
    electronic_integrator::E
end

OrdinaryDiffEq.isfsal(::BABwithTsit5) = false

# OrdinaryDiffEq.get_fsalfirstlast(cache::BABwithTsit5Cache, u::Any) = (nothing, nothing)

function OrdinaryDiffEq.alg_cache(
    alg::BABwithTsit5,u,rate_prototype,
    ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits},
    ::Type{tTypeNoUnits},
    uprev,
    uprev2,
    f,
    t,
    dt,
    reltol,
    p,
    calck,
    inplace::Val{true}
    ) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    
    tmp = zero(u)
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))

    electronic_problem = DynamicsMethods.DensityMatrixODEProblem(
        Array(DynamicsUtils.get_quantum_subsystem(u)),
        (0.0, dt),
        NQCModels.nstates(p)
    )
    electronic_integrator = SciMLBase.init(
        electronic_problem, alg.electronic_algorithm;
        save_on=false, save_everystep=false, dt=dt/5, adaptive=false
    )

    BABwithTsit5Cache(u, uprev, tmp, vtmp, k, electronic_integrator)
end

function OrdinaryDiffEq.initialize!(integrator, integrator_cache::BABwithTsit5Cache)

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

@muladd function OrdinaryDiffEq.perform_step!(integrator, cache::BABwithTsit5Cache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, vtmp, electronic_integrator = cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)
    σprev = DynamicsUtils.get_quantum_subsystem(uprev)

    rfinal = DynamicsUtils.get_positions(u)
    vfinal = DynamicsUtils.get_velocities(u)

    step_B!(vtmp, vprev, dt/2, k)
    step_A!(rfinal, rprev, dt, vtmp)

    NQCCalculators.update_cache!(p.cache, rfinal)
    if p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
       DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t, σprev)
    elseif p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t, p.method.state)
    end
    step_B!(DynamicsUtils.get_velocities(u), vtmp, dt/2, k)

    d = NQCCalculators.get_nonadiabatic_coupling(p.cache, rfinal)
    vals = NQCCalculators.get_eigen(p.cache, rfinal).values
    DynamicsMethods.update_parameters!(electronic_integrator.p, vals, d, vfinal, t+dt)

    density_matrix = DynamicsUtils.get_quantum_subsystem(u)
    set_ut!(electronic_integrator, density_matrix, t)
    SciMLBase.step!(electronic_integrator, dt, true)
    copy!(density_matrix, electronic_integrator.u)
end
