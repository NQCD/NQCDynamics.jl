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
        u .+= b
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
    HDF5.h5open(reduction.filename, "w") do file
        for (i, trajectory_id) in enumerate(I)
            trajectory_group = HDF5.create_group(file, "trajectory_$(trajectory_id)")
            trajectory = batch[i]
            for k in keys(trajectory)
                value = trajectory[k]
                if (value isa Array{<:Number}) || (value isa Number)
                    trajectory_group[string(k)] = trajectory[k]
                elseif (value isa Vector{<:Array})
                    output = cat(value...; dims=3)
                    trajectory_group[string(k)] = output
                else
                    throw(error("Cannot convert output type to HDF5 format"))
                end
            end
        end
    end
    return ("Output written to $(reduction.filename).", false)
end
