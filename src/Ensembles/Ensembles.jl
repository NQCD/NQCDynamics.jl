module Ensembles

using ..NonadiabaticMolecularDynamics
using ..Dynamics: AbstractMDEF
using ..InitialConditions: DynamicalDistribution
using SciMLBase: remake, EnsembleThreads, EnsembleProblem, solve
using DocStringExtensions
using RecursiveArrayTools: ArrayPartition
using TypedTables

function select_u0(::Simulation{<:Union{Classical, AbstractMDEF}}, v, r)
    ArrayPartition(v, r)
end

function select_u0(::RingPolymerSimulation{<:Classical}, v, r)
    RingPolymerClassicalDynamicals(v, r)
end

function select_u0(sim::Simulation{<:FSSH}, v, r)
    SurfaceHoppingVariables(v, r, sim.calculator.model.n_states, distribution.state)
end

function select_u0(sim::RingPolymerSimulation{<:NRPMD}, v, r)
    RingPolymerMappingVariables(v, r, sim.calculator.model.n_states, distribution.state)
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
    kwargs...
    )

    stripped_kwargs = austrip_kwargs(;kwargs...)

    problem = Dynamics.create_problem(select_u0(sim, selection.distribution, 1), austrip.(tspan), sim)
    problem = remake(problem, callback=Dynamics.get_callbacks(sim))
    ensemble_problem = EnsembleProblem(
        problem,
        prob_func=selection,
        output_func=output,
        reduction=reduction
    )

    solve(ensemble_problem, Dynamics.select_algorithm(sim), ensemble_algorithm; stripped_kwargs...)
end

function run_ensemble_standard_output(
    sim::AbstractSimulation, tspan,
    selection; output=(:u), ensemble_algorithm=EnsembleThreads(), kwargs...)

    stripped_kwargs = austrip_kwargs(;kwargs...)

    problem = Dynamics.create_problem(select_u0(sim, selection.distribution, 1), austrip.(tspan), sim)
    problem = remake(problem, callback=Dynamics.get_callbacks(sim))

    new_selection = SelectWithCallbacks(selection, Dynamics.get_callbacks(sim), output, kwargs[:trajectories])

    ensemble_problem = EnsembleProblem(problem, prob_func=new_selection)

    solve(ensemble_problem, Dynamics.select_algorithm(sim), ensemble_algorithm; stripped_kwargs...)
    [Table(t=vals.t, vals.saveval) for vals in new_selection.values]
end

end # module