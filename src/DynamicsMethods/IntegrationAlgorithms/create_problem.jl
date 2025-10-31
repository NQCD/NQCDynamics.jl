using StochasticDiffEq: BAOABCache, BAOABConstantCache
using SciMLBase
using NQCDynamics.DynamicsMethods.ClassicalMethods

function DynamicsMethods.create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF})
    StochasticDiffEq.DynamicalSDEProblem(ClassicalMethods.acceleration!, DynamicsUtils.velocity!, SciMLBase.DynamicalNoiseFunction(matrix_friction_update!),
        DynamicsUtils.get_velocities(u0), DynamicsUtils.get_positions(u0), tspan, sim)
end

function matrix_friction_update!(dutmp, utmp, integrator, cache::BAOABCache)
    @unpack t,dt,p = integrator
    @unpack half, gtmp = cache

    DynamicsMethods.ClassicalMethods.friction!(gtmp, utmp, p, t+dt*half)
    DynamicsMethods.IntegrationAlgorithms.step_O!(cache, integrator)
end

function matrix_friction_update!(utmp, integrator, cache::BAOABConstantCache)
    @unpack t,dt,sqdt,p,W = integrator
    @unpack half = cache
    Λ = integrator.g(utmp,p,t+dt*half) # friction tensor
    # noise strength: σ = square root of (temperature / mass) for each atom 
    σ = repeat(@. sqrt(get_temperature(p, t+dt*half) / p.atoms.masses);inner=ndofs(p))
    # eigen decomposition of Λ
    γ, c = LAPACK.syev!('V', 'U', Λ) # symmetric eigen (so far just calculated some values)
    clamp!(γ, 0, Inf) # makes positive semidefinite
    c1 = diagm(exp.(-γ.*dt))
    c2 = diagm(sqrt.(1 .- diag(c1).^2)) #more calculating shit 
    noise = σ.*W.dW[:] / sqdt # W.dW[:] noise rescaled by square root of time step to make it normal distributed
    du3 = c*c1*c'*du2[:] + c*c2*c'*noise # O step from Leimkuhler and Matthews (2013) 
    return reshape(du3, size(du2))
end