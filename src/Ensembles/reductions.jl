using StatsBase: mean
using HDF5: HDF5
using Dictionaries: Dictionary
using NQCBase: NQCBase

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
    XYZFileReduction(filename, sim)

Write trajectory positions to (ext)XYZ files, one file per trajectory.

Each trajectory is written to a separate file with the trajectory ID appended
to the base filename (e.g., `output_1.xyz`, `output_2.xyz`).

The `OutputPosition` output must be included in the `output` tuple passed to
[`run_dynamics`](@ref). If `OutputPosition` is not present, no file will be written
for that trajectory.

For simulations with a `PeriodicCell`, the extended XYZ format is used, including
cell and periodic boundary information. For other cell types (e.g., `InfiniteCell`),
a plain XYZ file is written with positions converted to Ångströms.

# Example
```julia
sim = Simulation(Atoms([:H]), Harmonic())
run_dynamics(sim, (0.0, 1.0), u0; output=OutputPosition, reduction=XYZFileReduction("output.xyz", sim))
# Writes output_1.xyz
```
"""
struct XYZFileReduction
    filename::String
    atoms
    cell

    function XYZFileReduction(filename::String, sim::AbstractSimulation)
        splitname = splitext(filename)
        ext = splitname[2]
        if ext in (".xyz", ".extxyz")
            return new(filename, sim.atoms, sim.cell)
        else
            return new(splitname[1] * ".xyz", sim.atoms, sim.cell)
        end
    end
end

function (reduction::XYZFileReduction)(u, batch, I)
    base, ext = splitext(reduction.filename)
    for (i, trajectory_id) in enumerate(I)
        trajectory = batch[i]
        if haskey(trajectory, :OutputPosition)
            positions = trajectory[:OutputPosition]
            filename = base * "_$(trajectory_id)" * ext
            _write_xyz_trajectory(filename, reduction.atoms, positions, reduction.cell)
        end
    end
    return ("Output written to $(base)_*$(ext).", false)
end

function _write_xyz_trajectory(filename, atoms, positions, cell::NQCBase.PeriodicCell)
    NQCBase.write_extxyz(filename, atoms, positions, cell)
end

function _write_xyz_trajectory(filename, atoms, positions, cell)
    open(filename, "w") do io
        for R in positions
            println(io, length(atoms))
            println(io, "")
            for (j, symbol) in enumerate(atoms.types)
                coords = NQCBase.au_to_ang.(R[:, j])
                println(io, "$(symbol) $(join(coords, " "))")
            end
        end
    end
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
