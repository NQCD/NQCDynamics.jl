using OrdinaryDiffEq.OrdinaryDiffEqCore: SciMLBase.AbstractDEAlgorithm, ODEIntegrator, loopfooter!, loopheader!, check_error!, perform_step!, postamble!
using OrdinaryDiffEq.OrdinaryDiffEqCore
using SciMLBase
using DiffEqBase
using SciMLBase: AbstractDEProblem, AbstractODEProblem, DEIntegrator, __init,  solve!

mutable struct CoupledODEProblem{P1<:SciMLBase.AbstractODEProblem, P2<:SciMLBase.AbstractODEProblem} <: SciMLBase.AbstractDEProblem
    prob1::P1
    prob2::P2
    callback::CallbackSet
end

mutable struct CoupledODEIntegrator{I1<:DiffEqBase.AbstractODEIntegrator, I2<:DiffEqBase.AbstractODEIntegrator, C}
    int1::I1
    int2::I2
    callback::C
end

struct CoupledODESolution{S1,S2}
    sol1::S1
    sol2::S2
end


# Make this as thin as possible - just init and dispatch
@inline function SciMLBase.__solve(
        prob::CoupledODEProblem,
        alg1::SciMLBase.AbstractDEAlgorithm, alg2::SciMLBase.AbstractDEAlgorithm, dt1::Real, dt2::Real, callback)
    
    integrator::CoupledODEIntegrator = SciMLBase.__init(prob, alg1, alg2, dt1, dt2, callback)
    
    return solve!(integrator)
end

function SciMLBase.__init(prob::CoupledODEProblem{P1,P2},
                          alg1::A1, alg2::A2,
                          dt1::Real, dt2::Real, callback) where
                          {P1<:SciMLBase.AbstractODEProblem, P2<:SciMLBase.AbstractODEProblem,
                           A1<:SciMLBase.AbstractDEAlgorithm, A2<:SciMLBase.AbstractDEAlgorithm} 
    
    int1 = SciMLBase.__init(prob.prob1, alg1; dt=dt1)
    int2 = SciMLBase.__init(prob.prob2, alg2; dt=dt2)
 
    return _init_coupled(int1, int2, callback)
end

@noinline function _init_coupled(int1::I1, int2::I2, callback::C) where {I1, I2, C}
    return CoupledODEIntegrator{I1,I2, C}(int1, int2, callback)
end

function SciMLBase.solve!(integrator::CoupledODEIntegrator{I1,I2}) where {I1, I2}
    @inbounds while !isempty(integrator.int1.opts.tstops) && !isempty(integrator.int2.opts.tstops)
        while integrator.int1.tdir * integrator.int1.t < first(integrator.int1.opts.tstops) && integrator.int2.tdir * integrator.int2.t < first(integrator.int2.opts.tstops)
            if (integrator.int1.do_error_check && check_error!(integrator.int1) != ReturnCode.Success) && (integrator.int2.do_error_check && check_error!(integrator.int2) != ReturnCode.Success)
                return CoupledODESolution(integrator.int1.sol, integrator.int2.sol)
            end
            OrdinaryDiffEqCore.perform_step!(integrator)
            if isempty(integrator.int1.opts.tstops) && isempty(integrator.int2.opts.tstops)
                break
            end
        end
        OrdinaryDiffEqCore.handle_tstop!(integrator.int1)
        OrdinaryDiffEqCore.handle_tstop!(integrator.int1)
    end
    postamble!(integrator.int1)
    postamble!(integrator.int2)
    if (integrator.int1.sol.retcode == ReturnCode.Default)
        integrator.int1.sol = SciMLBase.solution_new_retcode(integrator.int1.sol, ReturnCode.Success)
    end
    if (integrator.int2.sol.retcode == ReturnCode.Default)
        integrator.int2.sol = SciMLBase.solution_new_retcode(integrator.int2.sol, ReturnCode.Success)
    end
    return CoupledODESolution(integrator.int1.sol, integrator.int2.sol)
end

function OrdinaryDiffEqCore.perform_step!(integrator::CoupledODEIntegrator)
    dt = max(integrator.int1.dt, integrator.int2.dt)
    OrdinaryDiffEqCore.step!(integrator.int1, dt, true)
    OrdinaryDiffEqCore.step!(integrator.int2, dt, true)
    OrdinaryDiffEqCore._loopfooter!(integrator)
    integrator.int1.p[2] .= integrator.int2.u
    integrator.int2.p[2] .= integrator.int1.u
end

function DiffEqBase.solve(
        prob::CoupledODEProblem,
        alg1::SciMLBase.AbstractDEAlgorithm,
        alg2::SciMLBase.AbstractDEAlgorithm,
        dt1::Real, dt2::Real, 
        callback::CallbackSet=CallbackSet())

    return SciMLBase.__solve(prob, alg1, alg2, dt1, dt2, callback)
end

function OrdinaryDiffEqCore._loopfooter!(integrator::CoupledODEIntegrator)
    OrdinaryDiffEqCore.handle_callbacks!(integrator)
end

function OrdinaryDiffEqCore.handle_callbacks!(integrator::CoupledODEIntegrator)
    continous_callback = integrator.callback.continuous_callbacks
    discrete_callback = integrator.callback.discrete_callbacks
    if !(continous_callback isa Tuple{})
        for i in eachindex(continous_callback)
            t = integrator.int1.t
            if continous_callback[i].condition(integrator.int1.u, integrator.int2.u, t, integrator)
                continous_callback[i].affect!(integrator.int1.u, integrator.int2.u, integrator)
            end
        end
    end
    if !(discrete_callback isa Tuple{})
        for i in eachindex(discrete_callback)
            t = integrator.int1.t
            if discrete_callback[i].condition(integrator.int1.u, integrator.int2.u, t, integrator)
                discrete_callback[i].affect!(integrator.int1.u, integrator.int2.u, integrator)
            end
        end
    end
end

function DiffEqBase.solve(prob::SciMLBase.EnsembleProblem{<:CoupledODEProblem}, args...; kwargs...)
    alg = extract_coupled_alg(args, kwargs, kwargs)

    if length(args) > 1
        SciMLBase.__solve(prob, alg, Base.tail(args)...; kwargs...)
    else
        SciMLBase.__solve(prob, alg; kwargs...)
    end
end

@inline function extract_coupled_alg(solve_args, solve_kwargs, prob_kwargs)
    if isempty(solve_args) || isnothing(first(solve_args))
        if haskey(solve_kwargs, :alg)
            solve_kwargs[:alg]
        elseif haskey(prob_kwargs, :alg)
            prob_kwargs[:alg]
        else
            nothing
        end
    elseif first(solve_args) isa Tuple{SciMLBase.AbstractSciMLAlgorithm, SciMLBase.AbstractSciMLAlgorithm} &&
           !(first(solve_args) isa SciMLBase.EnsembleAlgorithm)
        first(solve_args)
    else
        nothing
    end
end

function DiffEqBase.solve(
        prob::CoupledODEProblem,
        alg;
        kwargs...)
    alg1, alg2 = alg
    dt1, dt2 = kwargs[:dt]
    callback = prob.callback
    return SciMLBase.__solve(prob, alg1, alg2, dt1, dt2, callback)
end