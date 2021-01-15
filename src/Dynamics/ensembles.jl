export run_ensemble
using ..InitialConditions

function run_ensemble(distribution::PhasespaceDistribution, tspan::Tuple, sim::AbstractSimulation; kwargs...)
    problem = create_problem(Phasespace(rand(distribution)...), tspan, sim)
    prob_func = get_problem_function(distribution)

    ensemble_problem = EnsembleProblem(problem, prob_func=prob_func)
    solve(ensemble_problem, select_algorithm(sim); kwargs...)
end

function get_problem_function(distribution::PhasespaceDistribution)
    prob_func(prob, i, repeat) = remake(prob, u0=Phasespace(rand(distribution)...))
end
