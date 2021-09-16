module Ensembles

using ..NonadiabaticMolecularDynamics: AbstractSimulation, Simulation, RingPolymerSimulation
using ..Dynamics:
    Dynamics,
    ClassicalMethods,
    SurfaceHoppingMethods,
    EhrenfestMethods,
    MappingVariableMethods
using ..InitialConditions: DynamicalDistribution

using SciMLBase: remake, EnsembleThreads, EnsembleProblem, solve
using DocStringExtensions
using RecursiveArrayTools: ArrayPartition
using TypedTables: TypedTables

function select_u0(sim::Simulation{<:Union{ClassicalMethods.Classical, ClassicalMethods.AbstractMDEF}}, v, r, state, type)
    Dynamics.DynamicsVariables(sim, v, r)
end

function select_u0(sim::RingPolymerSimulation{<:ClassicalMethods.ThermalLangevin}, v, r, state, type)
    Dynamics.DynamicsVariables(sim, RingPolymerArray(v), RingPolymerArray(r))
end

function select_u0(sim::AbstractSimulation{<:SurfaceHoppingMethods.FSSH}, v, r, state, type)
    return Dynamics.DynamicsVariables(sim, v, r, state; type=type)
end

function select_u0(sim::AbstractSimulation{<:EhrenfestMethods.Ehrenfest}, v, r, state, type)
    return Dynamics.DynamicsVariables(sim, v, r, state; type=type)
end

function select_u0(sim::RingPolymerSimulation{<:MappingVariableMethods.NRPMD}, v, r, state, type)
    Dynamics.DynamicsVariables(sim, v, r, state; type=type)
end

function select_u0(sim::AbstractSimulation{<:SurfaceHoppingMethods.IESH}, v, r, state, type)
    Dynamics.DynamicsVariables(sim, v, r)
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
