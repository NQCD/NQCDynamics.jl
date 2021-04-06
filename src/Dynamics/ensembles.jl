export run_ensemble
using ..InitialConditions

function run_ensemble(distribution::DynamicalDistribution,
                      tspan::Tuple, sim::AbstractSimulation;
                      output_func=(sol,i)->(sol,false),
                      reduction = (u,data,I)->(append!(u,data),false),
                      state::Integer=1,
                      output=(:u,),
                      callback=nothing,
                      kwargs...)

    stripped_kwargs = austrip_kwargs(;kwargs...)
    callbacks, values = create_saving_callbacks(kwargs[:trajectories], output)

    problem = create_problem(select_u0(sim, distribution, state), austrip.(tspan), sim)
    prob_func = get_problem_function(sim, distribution, state, callbacks, callback)

    ensemble_problem = EnsembleProblem(problem,
                                       prob_func=prob_func,
                                       reduction=reduction,
                                       output_func=output_func)
    solve(ensemble_problem, select_algorithm(sim); stripped_kwargs...)
    [Table(t=vals.t, vals.saveval) for vals in values]
end

function create_saving_callbacks(trajectories, output)
    callbacks = []
    values = []
    for i=1:trajectories
        cb, vals = create_saving_callback(output)
        push!(callbacks, cb)
        push!(values, vals)
    end
    callbacks, values
end

function get_problem_function(sim::AbstractSimulation, distribution::DynamicalDistribution, state::Integer, callbacks, callback)
    function prob_func(prob, i, repeat)
        callback_set = CallbackSet(callback, callbacks[i], get_callbacks(sim))
        remake(prob, u0=select_u0(sim, distribution, state), callback=callback_set)
    end
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
