using OrdinaryDiffEq.OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, ODEIntegrator, loopfooter!, loopheader!, check_error!, perform_step!, postamble!
using OrdinaryDiffEq.OrdinaryDiffEqCore
using SciMLBase
using DiffEqBase
using SciMLBase: AbstractDEProblem, AbstractODEProblem, DEIntegrator, __init,  solve!

mutable struct CoupledODEProblem <: SciMLBase.AbstractDEProblem
    prob1::SciMLBase.AbstractODEProblem
    prob2::SciMLBase.AbstractODEProblem
end

mutable struct CoupledODEIntegrator
    int1::ODEIntegrator
    int2::ODEIntegrator
end

function SciMLBase.__solve(
        prob::CoupledODEProblem,
        alg1::OrdinaryDiffEqAlgorithm, alg2::OrdinaryDiffEqAlgorithm, dt1::Real, dt2::Real, args...;
        kwargs...)
    integrator = SciMLBase.__init(prob, alg1, alg2, dt1, dt2, args...; kwargs...)
    solve!(integrator)
    (integrator.int1.sol, integrator.int2.sol)
end

function SciMLBase.__init(prob::CoupledODEProblem, alg1::OrdinaryDiffEqAlgorithm, alg2::OrdinaryDiffEqAlgorithm,
                          dt1::Real, dt2::Real, args...; kwargs...)
    int1 = SciMLBase.__init(prob.prob1, alg1, dt=dt1, args...; kwargs...)
    int2 = SciMLBase.__init(prob.prob2, alg2, dt=dt2, args...; kwargs...)
    integrator = CoupledODEIntegrator(int1, int2)
end

function SciMLBase.solve!(integrator::CoupledODEIntegrator)
    @inbounds while !isempty(integrator.int1.opts.tstops) || !isempty(integrator.int2.opts.tstops)
        while integrator.int1.tdir * integrator.int1.t < first(integrator.int1.opts.tstops) || integrator.int2.tdir * integrator.int2.t < first(integrator.int2.opts.tstops)
            if (integrator.int1.do_error_check && check_error!(integrator.int1) != ReturnCode.Success) && (integrator.int2.do_error_check && check_error!(integrator.int2) != ReturnCode.Success)
                return (integrator.int1.sol, integrator.int2.sol)
            end
            OrdinaryDiffEqCore.perform_step!(integrator)
            if isempty(integrator.int1.opts.tstops) && isempty(integrator.int2.opts.tstops)
                break
            end
        end
    end

    postamble!(integrator.int1)
    postamble!(integrator.int2)

    if (integrator.int1.sol.retcode != ReturnCode.Default) && (integrator.int2.sol.retcode != ReturnCode.Default)
        return (integrator.int1.sol, integrator.int2.sol)
    end
    (integrator.int1.sol, integrator.int2.sol) = (SciMLBase.solution_new_retcode(integrator.int1.sol, ReturnCode.Success), SciMLBase.solution_new_retcode(integrator.int2.sol, ReturnCode.Success))
end

function OrdinaryDiffEqCore.perform_step!(integrator::CoupledODEIntegrator)
    dt = max(integrator.int1.dt, integrator.int2.dt)
    OrdinaryDiffEqCore.step!(integrator.int1, dt, true)
    OrdinaryDiffEqCore.step!(integrator.int2, dt, true)
    integrator.int1.p[2] .= integrator.int2.cache.u
    integrator.int2.p[2] .= integrator.int1.cache.u
end

function DiffEqBase.solve(
        prob::CoupledODEProblem,
        args...;
        alg1::OrdinaryDiffEqAlgorithm, alg2::OrdinaryDiffEqAlgorithm, dt1::Real, dt2::Real,
        kwargs...)

    SciMLBase.__solve(prob, alg1, alg2, dt1, dt2, args...; kwargs...)
end