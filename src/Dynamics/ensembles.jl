export run_ensemble
export run_dissociation_ensemble
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

function select_u0(sim::Simulation{<:IESH}, distribution::DynamicalDistribution, state::Vector{Int})
    SurfaceHoppingDynamicals(rand(distribution)..., sim.calculator.model.n_states, state)
end

function select_u0(sim::RingPolymerSimulation{<:NRPMD}, distribution::DynamicalDistribution, state::Integer)
    RingPolymerMappingDynamicals(rand(distribution)..., sim.calculator.model.n_states, state)
end

function run_dissociation_ensemble(
    distribution::DynamicalDistribution,
    tspan::Tuple, sim::AbstractSimulation;
    atom_indices=(1,2),
    distance=1.0,
    kwargs...)

    stripped_kwargs = austrip_kwargs(;kwargs...)

    u0 = ArrayPartition(InitialConditions.pick(distribution, 1)...)
    problem = create_problem(u0, austrip.(tspan), sim)

    @views function atoms_far_apart(u, t, integrator)
        R = get_positions(u)
        norm(R[:,atom_indices[1]] .- R[:,atom_indices[2]]) > austrip(distance)
    end
    terminate = create_terminating_callback(atoms_far_apart)

    function output_func(sol, i)
        if atoms_far_apart(sol.u[end], 0.0, 0.0)
            return (1, false)
        else
            return (0, false)
        end
    end

    function prob_func(prob, i, repeat)
        u0 = ArrayPartition(InitialConditions.pick(distribution, i)...)
        remake(prob, u0=u0)
    end

    reduction(u,batch,I) = (mean(batch), false)

    ensemble_problem = EnsembleProblem(problem,
                                       prob_func=prob_func,
                                       reduction=reduction,
                                       output_func=output_func)

    solve(ensemble_problem, select_algorithm(sim); callback=terminate, stripped_kwargs...)
end
