using RingPolymerArrays: RingPolymerArrays

struct BCBwithTsit5 <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

OrdinaryDiffEq.isfsal(::BCBwithTsit5) = false

mutable struct BCBwithTsit5Cache{uType,rType,vType,rateType,uEltypeNoUnits,T} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    rtmp::rType
    vtmp::vType
    k::rateType
    cayley::Vector{Matrix{uEltypeNoUnits}}
    tsit5cache::T
end

function OrdinaryDiffEq.alg_cache(::BCBwithTsit5,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,inplace::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    σ = zero(DynamicsUtils.get_quantum_subsystem(u))
    σrate = zero(DynamicsUtils.get_quantum_subsystem(rate_prototype))
    σprev = zero(DynamicsUtils.get_quantum_subsystem(uprev))
    σprev2 = zero(DynamicsUtils.get_quantum_subsystem(uprev2))
    tsit5cache = OrdinaryDiffEq.alg_cache(OrdinaryDiffEq.Tsit5(), σ, σrate, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, σprev, σprev2, f, t, dt, reltol, p, calck, inplace)
    tmp = zero(u)
    rtmp = zero(DynamicsUtils.get_positions(u))
    vtmp = zero(DynamicsUtils.get_velocities(u))
    k = zero(DynamicsUtils.get_positions(rate_prototype))
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=false)
    BCBwithTsit5Cache(u, uprev, tmp, rtmp, vtmp, k, cayley, tsit5cache)
end

function OrdinaryDiffEq.initialize!(integrator, cache::BCBwithTsit5Cache)

    r = DynamicsUtils.get_positions(integrator.u)
    v = DynamicsUtils.get_velocities(integrator.u)
    σprev = DynamicsUtils.get_quantum_subsystem(integrator.u)
    Calculators.update_electronics!(integrator.p.calculator, r)
    if integrator.p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsMethods.EhrenfestMethods.acceleration!(cache.k, v, r, integrator.p, integrator.t, σprev)
    elseif integrator.p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsMethods.SurfaceHoppingMethods.acceleration!(cache.k, v, r, integrator.p, integrator.t, integrator.p.method.state)
    end
end

@muladd function OrdinaryDiffEq.perform_step!(integrator, cache::BCBwithTsit5Cache, repeat_step=false)
    @unpack t, dt, uprev, u, p = integrator
    @unpack k, rtmp, vtmp, cayley = cache

    rprev = DynamicsUtils.get_positions(uprev)
    vprev = DynamicsUtils.get_velocities(uprev)
    σprev = DynamicsUtils.get_quantum_subsystem(uprev)

    copyto!(rtmp, rprev)

    step_B!(vtmp, vprev, dt/2, k)
    RingPolymerArrays.transform_to_normal_modes!(rtmp, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(vtmp, p.beads.transformation)
    step_C!(vtmp, rtmp, cayley)
    RingPolymerArrays.transform_from_normal_modes!(rtmp, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(vtmp, p.beads.transformation)

    Calculators.update_electronics!(p.calculator, rtmp)
    if p.method isa DynamicsMethods.EhrenfestMethods.AbstractEhrenfest
        DynamicsMethods.EhrenfestMethods.acceleration!(k, vtmp, rtmp, p, t, σprev)
    elseif p.method isa DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping
        DynamicsMethods.SurfaceHoppingMethods.acceleration!(k, vtmp, rtmp, p, t, p.method.state)
    end
    step_B!(DynamicsUtils.get_velocities(u), vtmp, dt/2, k)

    copyto!(DynamicsUtils.get_positions(u), rtmp)

    # Tsit5 algorithm modified for only quantum subsystem
    @unpack a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,a62,a63,a64,a65,a71,a72,a73,a74,a75,a76 = cache.tsit5cache.tab
    @unpack k1,k2,k3,k4,k5,k6,k7,tmp = cache.tsit5cache

    DynamicsUtils.set_quantum_derivative!(k1, vtmp, σprev, p)
    a = dt*a21
    @.. tmp = σprev+a*k1
    DynamicsUtils.set_quantum_derivative!(k2, vtmp, tmp, p)
    @.. tmp = σprev+dt*(a31*k1+a32*k2)
    DynamicsUtils.set_quantum_derivative!(k3, vtmp, tmp, p)
    @.. tmp = σprev+dt*(a41*k1+a42*k2+a43*k3)
    DynamicsUtils.set_quantum_derivative!(k4, vtmp, tmp, p)
    @.. tmp = σprev+dt*(a51*k1+a52*k2+a53*k3+a54*k4)
    DynamicsUtils.set_quantum_derivative!(k5, vtmp, tmp, p)
    @.. tmp = σprev+dt*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)
    DynamicsUtils.set_quantum_derivative!(k6, vtmp, tmp, p)
    @.. tmp = σprev+dt*(a71*k1+a72*k2+a73*k3+a74*k4+a75*k5+a76*k6)
    copyto!(DynamicsUtils.get_quantum_subsystem(u), tmp)

    integrator.destats.nf += 6

end

DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping}) = BCBwithTsit5()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.EhrenfestMethods.AbstractEhrenfest}) = BCBwithTsit5()
