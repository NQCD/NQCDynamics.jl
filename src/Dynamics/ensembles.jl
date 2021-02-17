export run_ensemble
using ..InitialConditions

function run_ensemble(distribution::DynamicalDistribution,
                      tspan::Tuple, sim::AbstractSimulation;
                      output_func=(sol,i)->(sol,false),
                      reduction = (u,data,I)->(append!(u,data),false),
                      state::Integer=1,
                      kwargs...)
    problem = create_problem(select_u0(sim, distribution, state), tspan, sim)
    prob_func = get_problem_function(sim, distribution, state)

    ensemble_problem = EnsembleProblem(problem,
                                       prob_func=prob_func,
                                       reduction=reduction,
                                       output_func=output_func)
    solve(ensemble_problem, select_algorithm(sim); kwargs...)
end

function get_problem_function(sim::AbstractSimulation, distribution::DynamicalDistribution, state::Integer)
    prob_func(prob, i, repeat) = remake(prob, u0=select_u0(sim, distribution, state))
end

function select_u0(::Simulation{<:Union{Classical, AbstractMDEF}}, distribution::DynamicalDistribution, ::Integer)
    ArrayPartition(rand(distribution)...)
end

function select_u0(::RingPolymerSimulation{<:Classical}, distribution::DynamicalDistribution, ::Integer)
    RingPolymerClassicalDynamicals(rand(distribution)...)
end

function select_u0(sim::Simulation{<:FSSH}, distribution::DynamicalDistribution, state::Integer)
    SurfaceHoppingDynamicals(rand(distribution)..., sim.calculator.model.n_states, state)
end

function select_u0(sim::RingPolymerSimulation{<:NRPMD}, distribution::DynamicalDistribution, state::Integer)
    RingPolymerMappingDynamicals(rand(distribution)..., sim.calculator.model.n_states, state)
end
