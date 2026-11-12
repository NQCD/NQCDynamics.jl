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

struct SortByTrajectoryReduction end
(::SortByTrajectoryReduction)(u,data,I) = (append!(u,data), false)

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
                write_hdf5!(trajectory_group, string(k), trajectory[k])
            end
        end
    end
    return ("Output written to $(reduction.filename).", false)
end

"""
    write_hdf5!(parent, name, value)

Recursively write `value` into the HDF5 `parent` under `name`.
Numeric arrays and numbers are written directly, vectors of equally sized arrays are
stacked along a new last dimension, and dictionaries, structs and other vectors become
nested groups.
"""
function write_hdf5!(parent, name::AbstractString, value)
    if value isa Union{Number,AbstractString,Symbol}
        parent[name] = value isa Symbol ? string(value) : value
    elseif value isa AbstractArray{<:Number}
        parent[name] = collect(value)
    elseif value isa AbstractVector{<:AbstractArray{<:Number}} && !isempty(value) &&
           all(==(size(first(value))) ∘ size, value)
        parent[name] = reshape(reduce(hcat, vec.(value)), size(first(value))..., :)
    elseif value isa AbstractDict || value isa Dictionary
        group = HDF5.create_group(parent, name)
        for k in keys(value)
            write_hdf5!(group, string(k), value[k])
        end
    elseif value isa AbstractArray
        group = HDF5.create_group(parent, name)
        for (i, v) in enumerate(value)
            write_hdf5!(group, string(i), v)
        end
    elseif value isa Tuple
        group = HDF5.create_group(parent, name)
        for (i, v) in enumerate(value)
            write_hdf5!(group, string(i), v)
        end
    elseif isstructtype(typeof(value)) && fieldcount(typeof(value)) > 0
        group = HDF5.create_group(parent, name)
        for f in fieldnames(typeof(value))
            write_hdf5!(group, string(f), getfield(value, f))
        end
    else
        throw(error("Cannot convert output type $(typeof(value)) to HDF5 format"))
    end
    return nothing
end

"""
Organize outputs by output type rather than by trajectory.

Returns a `Dictionary` where keys are output names (e.g., `:OutputPosition`, `:OutputVelocity`)
and values are `Vector`s containing the corresponding output from each trajectory.

This makes it easier to compare a specific output across all trajectories without
needing to manually extract and collect data from individual trajectory dictionaries.

# Example
```julia
# With SortByTrajectoryReduction (default):
# results = [traj1_dict, traj2_dict, traj3_dict]
# traj1_dict = Dictionary(:Time => ..., :OutputPosition => ..., :OutputVelocity => ...)
# traj2_dict = Dictionary(:Time => ..., :OutputPosition => ..., :OutputVelocity => ...)
# traj3_dict = Dictionary(:Time => ..., :OutputPosition => ..., :OutputVelocity => ...)

# With SortByOutputReduction:
# results = Dictionary(
#     :OutputPosition => [traj1_pos, traj2_pos, traj3_pos],
#     :OutputVelocity => [traj1_vel, traj2_vel, traj3_vel],
#     :Time => [traj1_time, traj2_time, traj3_time]
# )
```
"""
mutable struct SortByOutputReduction
    initialised::Bool
    SortByOutputReduction() = new(false)
end

function (reduction::SortByOutputReduction)(u, batch, I)
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
