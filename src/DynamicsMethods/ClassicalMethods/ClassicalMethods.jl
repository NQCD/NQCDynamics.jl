
module ClassicalMethods

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsUtils,
    DynamicsMethods,
    Estimators,
    RingPolymers
using NQCCalculators
using NQCBase: Atoms
using NQCModels: Model, ndofs
using UnPack: @unpack
using FastBroadcast: @..

include("classical.jl")
export Classical
include("langevin.jl")
export Langevin, ThermalLangevin
include("mdef.jl")
export MDEF
export DiabaticMDEF
include("md_mdef.jl")
export MD_MDEF

include("rpmdef.jl")

function step_X_classical!(integrator, integrator_cache)
    return nothing
end

function step_X_constantfriction!(integrator, integrator_cache)
    @unpack t, dt, uprev, p, W, f = integrator
    @unpack tmp, half, gtmp, flatdutmp, tmp1, tmp2, noise, c1, c2 = integrator_cache

    r = DynamicsUtils.get_positions(tmp)
    v = DynamicsUtils.get_velocities(tmp)

    f.g(gtmp,r,p,t+dt*half) #friction function
    @. noise = gtmp .* W.dW[:]
    @.. v = c1 * v + c2 * noise
end

const ClassicalMethodUnion = Union{Classical, Langevin, ThermalLangevin, MDEF, DiabaticMDEF}

end # module
