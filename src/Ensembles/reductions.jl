using StatsBase: mean
using HDF5: HDF5
using Dictionaries: Dictionary

"""
Sum the outputs from each trajectory.
"""
mutable struct SumReduction
    initialised::Bool
    SumReduction() = new(false)
end
function (reduction::SumReduction)(u, batch, I)
    if !reduction.initialised
        u = get_initial_u(batch)
        reduction.initialised = true
    end
    sum_outputs!(u, batch)
    return (u, false)
end

function get_initial_u(batch)
    template = first(batch)
    return Dictionary(keys(template), [zero.(t) for t in template])
end

function sum_outputs!(u, batch)
    for b in batch
        mergewith!(+, u, b)
        # u .+= b
    end
end

"""
Average the outputs over all trajectories.
"""
mutable struct MeanReduction
    counter::Int
    initialised::Bool
    MeanReduction() = new(0, false)
end

function (reduction::MeanReduction)(u, batch, I)
    if !reduction.initialised
        u = get_initial_u(batch)
        reduction.initialised = true
    end
    u .*= reduction.counter
    reduction.counter += length(batch)
    sum_outputs!(u, batch)
    u ./= reduction.counter
    return (u, false)
end

struct AppendReduction end
(::AppendReduction)(u,data,I) = (append!(u,data), false)

struct FileReduction
    filename::String

    function FileReduction(filename::String)
        splitname = splitext(filename)
        ext = splitname[2]
        if (ext == ".h5") || (ext == ".hdf5")
            return new(filename)
        else
            filename = splitname[1] * ".h5"
            return new(filename)
        end
    end
end

function (reduction::FileReduction)(u, batch, I)
    HDF5.h5open(reduction.filename, "cw") do file
        for (i, trajectory_id) in enumerate(I)
            trajectory_group = HDF5.create_group(file, "trajectory_$(trajectory_id)")
            trajectory = batch[i]
            for k in keys(trajectory)
                value = trajectory[k]
                if (value isa Array{<:Number}) || (value isa Number)
                    trajectory_group[string(k)] = trajectory[k]
                elseif (value isa Vector{<:Array})
                    output = reshape(reduce(hcat, value), size(value[1])..., :)
                    trajectory_group[string(k)] = output
                else
                    throw(error("Cannot convert output type to HDF5 format"))
                end
            end
        end
    end
    return ("Output written to $(reduction.filename).", false)
end

"""
Organize outputs by output type rather than by trajectory.

Returns a `Dictionary` where keys are output names (e.g., `:OutputPosition`, `:OutputVelocity`)
and values are `Vector`s containing the corresponding output from each trajectory.

This makes it easier to compare a specific output across all trajectories without
needing to manually extract and collect data from individual trajectory dictionaries.

# Example
```julia
# With AppendReduction (default):
# results = [traj1_dict, traj2_dict, traj3_dict]
# traj1_dict = Dictionary(:Time => ..., :OutputPosition => ..., :OutputVelocity => ...)
# traj2_dict = Dictionary(:Time => ..., :OutputPosition => ..., :OutputVelocity => ...)
# traj3_dict = Dictionary(:Time => ..., :OutputPosition => ..., :OutputVelocity => ...)

# With OutputReduction:
# results = Dictionary(
#     :OutputPosition => [traj1_pos, traj2_pos, traj3_pos],
#     :OutputVelocity => [traj1_vel, traj2_vel, traj3_vel],
#     :Time => [traj1_time, traj2_time, traj3_time]
# )
```
"""
mutable struct OutputReduction
    initialised::Bool
    OutputReduction() = new(false)
end

function (reduction::OutputReduction)(u, batch, I)
    if !reduction.initialised
        u = get_initial_output_dict(batch)
        reduction.initialised = true
    end
    append_outputs!(u, batch)
    return (u, false)
end

function get_initial_output_dict(batch)
    template = first(batch)
    return Dictionary(keys(template), [typeof(v)[] for v in template])
end

function append_outputs!(u, batch)
    for trajectory in batch
        for k in keys(trajectory)
            push!(u[k], trajectory[k])
        end
    end
end
