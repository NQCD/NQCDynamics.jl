module Ensembles

using ..NonadiabaticMolecularDynamics: AbstractSimulation, Dynamics, austrip_kwargs
using ..InitialConditions: DynamicalDistribution
using SciMLBase: remake, EnsembleThreads, EnsembleProblem, solve
using DocStringExtensions

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

    problem = Dynamics.create_problem(Dynamics.select_u0(sim, selection.distribution, 1), austrip.(tspan), sim)
    ensemble_problem = EnsembleProblem(
        problem,
        prob_func=selection,
        output_func=output,
        reduction=reduction,
        # u_init=reduction.u_init
    )

    solve(ensemble_problem, Dynamics.select_algorithm(sim), ensemble_algorithm; stripped_kwargs...)
end


end # module