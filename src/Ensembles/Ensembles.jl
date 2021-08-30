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
    return DynamicsVariables(sim, v, r, state; type=type)
end

function select_u0(sim::AbstractSimulation{<:Ehrenfest}, v, r, state, type)
    return DynamicsVariables(sim, v, r, state; type=type)
end

function select_u0(sim::RingPolymerSimulation{<:NRPMD}, v, r, state, type)
    DynamicsVariables(sim, v, r, state; type=type)
end

function select_u0(sim::AbstractSimulation{<:IESH}, v, r, state, type)
    DynamicsVariables(sim, v, r)
end

include("selections.jl")
include("reductions.jl")
include("outputs.jl")

function run_ensemble(
    sim::AbstractSimulation, tspan,
    distribution;
    selection=nothing,
    output=(sol,i)->(sol,false),
    reduction=(u,data,I)->(append!(u,data),false),
    ensemble_algorithm=EnsembleThreads(),
    algorithm=Dynamics.select_algorithm(sim),
    kwargs...
    )

    stripped_kwargs = austrip_kwargs(;kwargs...)

    u0 = select_u0(sim, rand(distribution)..., distribution.state, distribution.type)
    problem = Dynamics.create_problem(u0, austrip.(tspan), sim)

    if hasfield(typeof(reduction), :u_init)
        u_init = reduction.u_init
    else
        u_init = []
    end

    if selection isa AbstractVector
        selection = OrderedSelection(distribution, selection)
    else
        selection = RandomSelection(distribution)
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

function run_ensemble_standard_output(sim::AbstractSimulation, tspan, distribution;
    selection=nothing, output=(:u),
    ensemble_algorithm=EnsembleThreads(), saveat=[], kwargs...)

    stripped_kwargs = austrip_kwargs(;kwargs...)
    saveat = austrip.(saveat)

    problem = Dynamics.create_problem(
        select_u0(sim, rand(distribution)...,
            distribution.state, distribution.type),
        austrip.(tspan),
        sim)

    if selection isa AbstractVector
        selection = OrderedSelection(distribution, selection)
    else
        selection = RandomSelection(distribution)
    end

    new_selection = SelectWithCallbacks(selection, Dynamics.get_callbacks(sim), output, kwargs[:trajectories], saveat=saveat)

    ensemble_problem = EnsembleProblem(problem, prob_func=new_selection)

    solve(ensemble_problem, Dynamics.select_algorithm(sim), ensemble_algorithm; saveat=saveat, stripped_kwargs...)
    [Table(t=vals.t, vals.saveval) for vals in new_selection.values]
end

end # module
