# debug function to print contents of a struct (particularly for NCQCalculators caches)
function print_contents(container)
    println("Contents of $(typeof(container)):\n")
    for name in fieldnames(typeof(container))
        println("$(name) = ", getfield(container, name))
        println("")
    end
end

"""
    run_dynamics(sim::AbstractSimulation, tspan, distribution;
        output,
        selection::Union{Nothing,AbstractVector}=nothing,
        reduction=AppendReduction(),
        ensemble_algorithm=SciMLBase.EnsembleSerial(),
        algorithm=DynamicsMethods.select_algorithm(sim),
        trajectories=1,
        kwargs...
        )

Run trajectories for timespan `tspan` sampling from `distribution`.

# Keywords

* `output` either a single function or a Tuple of functions with the signature `f(sol, i)` that takes the DifferentialEquations solution and returns the desired output quantity.
* `selection` should be an `AbstractVector` containing the indices to sample from the `distribution`. By default, `nothing` leads to random sampling.
* `reduction` defines how the data is reduced across trajectories. Options are `AppendReduction()`, `MeanReduction()`, `SumReduction` and `FileReduction(filename)`.
* `ensemble_algorithm` is the algorithm from DifferentialEquations which determines which form of parallelism is used.
* `algorithm` is the algorithm used to integrate the equations of motion.
* `trajectories` is the number of trajectories to perform.
* `savetime` saves the simulation time-steps, useful when usings with integrators with variable time-stepping.
* `precompile_dynamics` runs a single step of the dynamics simulation in order to optimise the full run, exploits a strange issue with the Julia compiler.
* `kwargs...` any additional keywords are passed to DifferentialEquations `solve``.
"""
function run_dynamics(
    sim::AbstractSimulation,
    tspan,
    distribution;
    output,
    selection::Union{Nothing,AbstractVector} = nothing,
    reduction = AppendReduction(),
    ensemble_algorithm = SciMLBase.EnsembleSerial(),
    algorithm = DynamicsMethods.select_algorithm(sim),
    trajectories = 1,
    savetime = true,
    precompile_dynamics = true,
    rjm = false,
    kwargs...,
)

    """
    Due to how Julia compiles code, running the same dynamics for a single time step forces precompilation, which will make subsequent simulations faster. 
    This effect is negligible for smaller models, but begins to matter for more computationally intensive dynamics. 
    
    See this GitHub issue for discussion: https://github.com/NQCD/NQCDynamics.jl/issues/365
    """
    function dynamics(tspan)
        @debug "Setting up dynamics with" tspan = tspan
        if !(output isa Tuple)
            output = (output,)
        end
        
        kwargs = NQCBase.austrip_kwargs(; kwargs...)
        trajectories = convert(Int, trajectories)
        tspan = austrip.(tspan)
        prob_func = Selection(distribution, selection, trajectories)
        
        u0 = sample_distribution(sim, distribution)
        @debug print_contents(sim.cache)
        problem = DynamicsMethods.create_problem(u0, tspan, sim)
        
        output_func = EnsembleSaver(output, savetime)
        
        ensemble_problem = SciMLBase.EnsembleProblem(
            problem;
            prob_func = prob_func,
            output_func = output_func,
            reduction = deepcopy(reduction), # deepcopy reduction function because it's used twice
        )
        
        return @timed SciMLBase.solve(
            ensemble_problem,
            algorithm,
            ensemble_algorithm;
            trajectories,
            kwargs...,
        )
    end

    if precompile_dynamics
        @debug "Beginning to precompile dynamics"
        # Get a short time limit that runs through saving at least one time step of data
        short_time = get( 
            kwargs, 
            :saveat, # Saveat will be an interval or a list. 
            get(
                kwargs,
                :dt,
                1.0 # dt should be a number / unitful qty
            )
        ) |> first # Get the number or the first list entry. 
        tspan_short = (0.0, short_time)
        precompile_time = dynamics(tspan_short)
        @debug "Reduction is" reduction=reduction
        if isa(reduction, FileReduction)
            @debug "FileReduction wrote files, which will now be deleted. "
            rm(reduction.filename)
        end
        @info "Pre-compiled dynamics in $(precompile_time.time) seconds."
    end

    if trajectories == 1
        @info "Performing 1 trajectory."
    else
        @info "Performing $trajectories trajectories."
    end

    stats = dynamics(tspan)
    log_simulation_duration(stats.time)

    if trajectories == 1
        if rjm NQCCalculators.get_maurergroup() end
        return stats.value.u[1]
    else
        if rjm NQCCalculators.get_maurergroup() end
        return stats.value.u
    end
end

function log_simulation_duration(duration_seconds)
    if duration_seconds < 60
        @info "Finished after $duration_seconds seconds."
    elseif duration_seconds < 3600
        duration_minutes = duration_seconds / 60
        @info "Finished after $duration_minutes minutes."
    else
        duration_hours = duration_seconds / 3600
        @info "Finished after $duration_hours hours."
    end
end
