module Ensembles

using ..NonadiabaticMolecularDynamics
using ..Dynamics: AbstractMDEF
using ..InitialConditions: DynamicalDistribution
using SciMLBase: remake, EnsembleThreads, EnsembleProblem, solve
using DocStringExtensions
using RecursiveArrayTools: ArrayPartition
using TypedTables

function select_u0(::Simulation{<:Union{Classical, AbstractMDEF}}, v, r, state, type)
    ArrayPartition(v, r)
end

function select_u0(::RingPolymerSimulation{<:ThermalLangevin}, v, r, state, type)
    ArrayPartition(RingPolymerArray(v), RingPolymerArray(r))
end

function select_u0(sim::AbstractSimulation{<:FSSH}, v, r, state, type)
    if type == :adiabatic
        return SurfaceHoppingVariables(v, r, sim.calculator.model.n_states, state)
    elseif type == :diabatic
        return SurfaceHoppingVariables(sim, v, r, state)
    else
        throw(ArgumentError("$type is not a recognised type argument.
            Try `:adiabatic` or `diabatic`."))
    end
end

function select_u0(sim::AbstractSimulation{<:Ehrenfest}, v, r, state, type)
    if type == :adiabatic
        return EhrenfestVariables(v, r, sim.calculator.model.n_states, state)
    elseif type == :diabatic
        return EhrenfestVariables(sim, v, r, state)
    else
        throw(ArgumentError("$type is not a recognised type argument.
            Try `:adiabatic` or `diabatic`."))
    end
end

function select_u0(sim::RingPolymerSimulation{<:NRPMD}, v, r, state, type)
    RingPolymerMappingVariables(v, r, sim.calculator.model.n_states, state)
end

function select_u0(sim::AbstractSimulation{<:IESH}, v, r, state, type)
    SurfaceHoppingVariablesIESH(v, r, sim.calculator.model.n_states, sim.method.n_electrons)
end

include("selections.jl")
include("reductions.jl")
include("outputs.jl")

function run_ensemble(
    sim::AbstractSimulation, tspan,
    selection;
    output=(sol,i)->(sol,false),
    reduction=(u,data,I)->(append!(u,data),false),
    ensemble_algorithm=EnsembleThreads(),
    algorithm=Dynamics.select_algorithm(sim),
    kwargs...
    )

    stripped_kwargs = austrip_kwargs(;kwargs...)

    u0 = select_u0(sim, rand(selection.distribution)..., selection.distribution.state, selection.distribution.type)
    problem = Dynamics.create_problem(u0, austrip.(tspan), sim)
    problem = remake(problem, callback=Dynamics.get_callbacks(sim))

    if hasfield(typeof(reduction), :u_init)
        u_init = reduction.u_init
    else
        u_init = []
    end

    ensemble_problem = EnsembleProblem(
        problem,
        prob_func=selection,
        output_func=output,
        reduction=reduction,
        u_init=u_init
    )

    solve(ensemble_problem, algorithm, ensemble_algorithm; stripped_kwargs...)
end

function run_ensemble_standard_output(
    sim::AbstractSimulation, tspan,
    selection; output=(:u), ensemble_algorithm=EnsembleThreads(), saveat=[], kwargs...)

    stripped_kwargs = austrip_kwargs(;kwargs...)
    saveat = austrip.(saveat)

    problem = Dynamics.create_problem(
        select_u0(sim, rand(selection.distribution)...,
            selection.distribution.state, selection.distribution.type),
        austrip.(tspan),
        sim)

    problem = remake(problem, callback=Dynamics.get_callbacks(sim))

    new_selection = SelectWithCallbacks(selection, Dynamics.get_callbacks(sim), output, kwargs[:trajectories], saveat=saveat)

    ensemble_problem = EnsembleProblem(problem, prob_func=new_selection)

    solve(ensemble_problem, Dynamics.select_algorithm(sim), ensemble_algorithm; saveat=saveat, stripped_kwargs...)
    [Table(t=vals.t, vals.saveval) for vals in new_selection.values]
end

end # module
