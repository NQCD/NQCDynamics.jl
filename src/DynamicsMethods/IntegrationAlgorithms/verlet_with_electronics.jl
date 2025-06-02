using OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache, update_coefficients!
using SciMLBase: SciMLBase, set_ut!
using NQCDynamics: DynamicsUtils
using .DynamicsUtils: acceleration!, get_positions, get_velocities, get_quantum_subsystem
using NQCCalculators

OrdinaryDiffEq.isfsal(::VerletwithElectronics) = false

mutable struct VerletwithElectronicsCache{uType,vType,rateType} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    vtmp::vType
    k::rateType
end

function OrdinaryDiffEq.alg_cache(::VerletwithElectronics,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    VerletwithElectronicsCache(u, uprev, tmp, vtmp, k)
end

function OrdinaryDiffEq.initialize!(integrator, integrator_cache::VerletwithElectronicsCache)
    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    NQCCalculators.update_cache!(integrator.p.cache, r)
    if integrator.p.method isa DynamicsMethods.SurfaceHoppingMethods.AbstractIESH
        DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, integrator.p.method.state)
    elseif integrator.p.method isa DynamicsMethods.EhrenfestMethods.EhrenfestNA
        ψ = DynamicsUtils.get_quantum_subsystem(integrator.u)
        DynamicsUtils.acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, ψ)
    end
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, integrator_cache::VerletwithElectronicsCache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, vtmp = integrator_cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)
    σprev = DynamicsUtils.get_quantum_subsystem(uprev)

    rfinal = DynamicsUtils.get_positions(u)
    vfinal = DynamicsUtils.get_velocities(u)
    σfinal = DynamicsUtils.get_quantum_subsystem(u)

    step_B!(vtmp, vprev, dt/2, k)
    step_A!(rfinal, rprev, dt, vtmp)

    NQCCalculators.update_cache!(p.cache, rfinal)
    if integrator.p.method isa DynamicsMethods.SurfaceHoppingMethods.AbstractIESH
        DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t, p.method.state)
    elseif integrator.p.method isa DynamicsMethods.EhrenfestMethods.EhrenfestNA
        DynamicsUtils.acceleration!(k, vtmp, rfinal, p, t, σprev)
    end

    step_B!(vfinal, vtmp, dt/2, k)

    DynamicsUtils.propagate_wavefunction!(σfinal, σprev, vfinal, rfinal, p, dt)

end

struct VerletwithElectronics2{T<:OrdinaryDiffEqAlgorithm,K<:Base.Pairs} <: OrdinaryDiffEqAlgorithm
    electronic_algorithm::T
    kwargs::K
end

VerletwithElectronics2(electronic_algorithm; kwargs...) = VerletwithElectronics2(electronic_algorithm, kwargs)

OrdinaryDiffEq.isfsal(::VerletwithElectronics2) = false

mutable struct VerletwithElectronics2Cache{uType,vType,rateType,E} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    vtmp::vType
    k::rateType
    electronic_integrator::E
end

function OrdinaryDiffEq.alg_cache(alg::VerletwithElectronics2,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))

    electronic_problem = DynamicsMethods.ElectronicODEProblem(Array(DynamicsUtils.get_quantum_subsystem(u)[:,1]), (0.0, dt), nstates(p))
    electronic_integrator = SciMLBase.init(electronic_problem, alg.electronic_algorithm; save_on=false, save_everystep=false, alg.kwargs...)

    VerletwithElectronics2Cache(u, uprev, tmp, vtmp, k, electronic_integrator)
end

function OrdinaryDiffEq.initialize!(integrator, integrator_cache::VerletwithElectronics2Cache)
    r = get_positions(integrator.u)
    v = get_velocities(integrator.u)
    NQCCalculators.update_cache!(integrator.p.cache, r)
    acceleration!(integrator_cache.k, v, r, integrator.p, integrator.t, integrator.p.method.state)
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, integrator_cache::VerletwithElectronics2Cache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, vtmp, electronic_integrator = integrator_cache
    sim = p

    rprev = get_positions(uprev)
    vprev = get_velocities(uprev)

    rfinal = get_positions(u)
    vfinal = get_velocities(u)

    # Velocity verlet steps
    step_B!(vtmp, vprev, dt/2, k)
    step_A!(rfinal, rprev, dt, vtmp)
    acceleration!(k, vtmp, rfinal, p, t, p.method.state)
    step_B!(vfinal, vtmp, dt/2, k)

    d = NQCCalculators.get_nonadiabatic_coupling(sim.cache, rfinal)
    vals = NQCCalculators.get_eigen(sim.cache, rfinal).values
    DynamicsMethods.update_parameters!(electronic_integrator.p, vals, d, vfinal, t+dt)

    for i in axes(get_quantum_subsystem(u), 2)
        c = view(get_quantum_subsystem(u), :, i)
        set_ut!(electronic_integrator, c, t)
        SciMLBase.step!(electronic_integrator, dt, true) # Step electronic integrator to current timestep
        get_quantum_subsystem(u)[:,i] .= electronic_integrator.u
    end

end
