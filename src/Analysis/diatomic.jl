"""
Analysis functions for surface chemistry of diatomic molecules.
"""
module Diatomic
using NQCDynamics: AbstractSimulation, Simulation, get_positions, Structure
using NQCBase
using Unitful, UnitfulAtomic
using LinearAlgebra
using Statistics

# Surface normal projection
surface_normal_height(x::AbstractVector, surface_normal::AbstractVector = [0, 0, 1]) =
    norm(x .* normalize(surface_normal))


# New conditions as of 2025-10-27
# H-H distance of a simulation frame
function H_H_distance(
    frame::AbstractArray, # MD frame in a DynamicsVariables type
    indices::Vector{Int},
    simulation::AbstractSimulation, # Simulation object informing PBC
)
    return Structure.pbc_distance(get_positions(frame), indices..., simulation.cell)
end

function H2_surface_distance(
    frame::AbstractArray, # MD frame of DynamicsVariables type
    diatomic_indices::Vector{Int}, # H atom indices
    toplayer_indices::Vector{Int}, # Top surface layer
    simulation::AbstractSimulation; # Simulation informing PBC and atom types
    surface_normal::AbstractVector = [0.0, 0.0, 1.0], # Surface normal vector used to define heights within the unit cell. 
)
    positions = get_positions(frame)
    slab_height = mean([
        surface_normal_height(x, surface_normal) for
        x in eachcol(positions[:, toplayer_indices])
    ])
    h2_height = surface_normal_height(
        Structure.pbc_center_of_mass(positions, diatomic_indices..., simulation),
        surface_normal,
    )
    return h2_height - slab_height
end


"""
    surface_distance_condition(x::AbstractArray, indices::Vector{Int}, simulation::AbstractSimulation; surface_distance_threshold=5.0*u"Å")

    Checks that the diatomic molecule is at least `surface_distance_threshold` away from the highest substrate atom in the simulation.
"""
function surface_distance_condition(
    x::AbstractArray, # DynamicsVariables
    indices::Vector{Int}, # Which indices are the diatomic species? All others are counted as substrate. 
    simulation::AbstractSimulation;
    surface_distance_threshold = austrip(5.0 * u"Å"), # Condition is true when the centre of mass of the diatomic is farther than this distance from the highest other atom. 
    surface_normal = Float64[0, 0, 1], # Unit direction orthogonal to the surface
)
    molecule_position = surface_normal_height(
        Structure.pbc_center_of_mass(x, indices..., simulation),
        surface_normal,
    ) # Height of the molecule in the surface normal direction

    # Get the height of all substrate atoms in surface normal direction. 
    substrate_heights = [
        surface_normal_height(get_positions(x)[:, substrate_atom_id], surface_normal)
        for substrate_atom_id in symdiff(1:length(simulation.atoms.masses), indices)
    ]

    # Ignore substrate above molecule in case PBC wrapping puts one above the diatomic
    highest_z = max(substrate_heights[substrate_heights.≤molecule_position]...)

    if molecule_position - highest_z ≥ surface_distance_threshold
        @debug "Surface distance condition evaluated true."
        return true
    else
        return false
    end
end

"""
    com_velocity_condition(x::AbstractArray, indices::Vector{Int}, simulation::AbstractSimulation; surface_normal::AbstractVector=[0,0,1])

    Evaluates true if the centre of mass velocity vector of the diatomic molecule points to the surface.
"""
function com_velocity_condition(
    x::AbstractArray,
    indices::Vector{Int},
    simulation::AbstractSimulation;
    surface_normal::AbstractVector = [0, 0, 1],
)
    if dot(
        Structure.velocity_center_of_mass(x, indices..., simulation),
        normalize(surface_normal),
    ) > 0
        return false
    else
        return true
    end
end

"""
    close_approach_condition(x::AbstractArray, indices::Vector{Int}, simulation::AbstractSimulation; threshold = 1.5u"Å")

    Evaluate true if the diatomic bond length is below `threshold`.
"""
function close_approach_condition(
    x::AbstractArray,
    indices::Vector{Int},
    simulation::AbstractSimulation;
    threshold = austrip(1.5u"Å"),
)
    if austrip(Structure.pbc_distance(x, indices..., simulation)) ≤ threshold
        return true
    else
        return false
    end
end

"""
    get_desorption_frame(
    trajectory::AbstractVector,
    diatomic_indices::Vector{Int},
    simulation::AbstractSimulation;
    surface_normal::Vector = [0, 0, 1],
    surface_distance_threshold = austrip(4.0 * u"Å"),
    fallback_distance_threshold = austrip(1.5u"Å"),
)

Determines the index in a trajectory where surface desorption begins.

This is evaluated using two conditions:

1. In the trajectory, the diatomic must be `surface_distance_threshold` or further away from the highest other atom. (In `surface_normal` direction).

2. The transition state to desorption is assumed to happen when the H-H-distance decreases below σ + μ for the last time in a trajectory. 

If the second condition is never reached, the whole trajectory is saved for troubleshooting. 
"""
function get_desorption_frame(
    trajectory::AbstractVector,
    diatomic_indices::Vector{Int},
    simulation::AbstractSimulation;
    surface_normal::Vector = [0, 0, 1],
    surface_distance_threshold = austrip(4.0 * u"Å"),
)
    # Calculate H2 surface distances and check if there was a desorption
    surface_distance = Float64[]
    sizehint!(surface_distance, length(trajectory))
    for idx in eachindex(trajectory)
        pos = get_positions(trajectory[idx])
        non_H_heights = [surface_normal_height(atom, surface_normal) for atom in eachcol(pos[:, symdiff(1:size(pos)[2], diatomic_indices)])]
        _, highest_H_index = findmax(non_H_heights)
        push!(
            surface_distance,
            H2_surface_distance(
                trajectory[idx],
                diatomic_indices,
                [highest_H_index],
                simulation;
                surface_normal = surface_normal,
            )
        )
    end
    
    desorption_frame = findfirst(surface_distance .≥ surface_distance_threshold)
    if isnothing(desorption_frame)
        @debug "No desorption found in trajectory (H-surface distance never above set threshold"
        return nothing
    end
    
    # Work backwards from desorption frame to determine H-H distance until it becomes larger than H2-surface distance
    H_H_distances = Float64[]
    desorption_frame = 0

    @views for idx in Iterators.reverse(trajectory)
        push!(
            H_H_distances, 
            H_H_distance(
                trajectory,
                diatomic_indices,
                simulation,
            ) |> austrip
        )
        h_h_distance_mean = mean(H_H_distances)
        h_h_distance_std = std(H_H_distances)
        desorption_frame += 1
        if (h_h_distance_mean + h_h_distance_std) < last(H_H_distances)
            @debug "Transition state observed with bond length $(last(h_h_distances)) > μ+σ $(h_h_distance_mean + h_h_distance_std) at index $(desorption_frame)"
            return length(trajectory) - desorption_frame + 1
        end
    end
    @warn "Check trajectory - Desorption observed, but H-surface distance never fell below H-H bond length. Returning index 1 to output full trajectory. "
    return 1
end

function get_desorption_angle(
    trajectory::AbstractVector,
    indices::Vector{Int},
    simulation::AbstractSimulation;
    surface_normal = [0, 0, 1],
    surface_distance_threshold = 5.0 * u"Å",
)
    # First determine where the reaction occurred on the surface.
    desorption_frame = get_desorption_frame(
        trajectory,
        indices,
        simulation;
        surface_distance_threshold = surface_distance_threshold,
        surface_normal = surface_normal,
    )
    if isnothing(desorption_frame)
        @debug "No desorption event detected in trajectory"
        return nothing
    end
    # Determine the average centre of mass velocity to decrease error due to vibration and rotation orthogonal to true translational component.
    com_velocities = zeros(
        eltype(trajectory[1]),
        length(surface_normal),
        length(trajectory) - desorption_frame,
    )
    for i = 1:length(trajectory)-desorption_frame
        com_velocities[:, i] .= Structure.velocity_center_of_mass(
            trajectory[i+desorption_frame],
            indices[1],
            indices[2],
            simulation,
        )
    end
    average_velocity = mean(com_velocities; dims = 2)
    # Now convert into an angle by arccos((a•b)/(|a|*|b|))
    return Structure.angle_between(vec(average_velocity), surface_normal)
end

# 6×6 Cartesian to internal coordinate transformation.
function transform_to_internal_coordinates(
    to_transform::Matrix,
    config::Matrix,
    index1::Int,
    index2::Int,
    sim::Simulation,
)
    U = transform_U(config, index1, index2, sim) # Generate transformation matrix
    return transpose(U) * to_transform * U
end

function transform_from_internal_coordinates(
    to_transform::Matrix,
    config::Matrix,
    index1::Int,
    index2::Int,
    sim::Simulation,
)
    U_inv = inv(transform_U(config, index1, index2, sim))
    return transpose(U_inv) * to_transform * U_inv
end

function transform_r(config, index1::Int, index2::Int)
    return sqrt(sum(mapslices(x -> (x[2] - x[1])^2, config[:, [index1, index2]]; dims = 2)))
end

function transform_r1(config, index1::Int, index2::Int)
    return transform_r(config[1:2, :], index1, index2)
end


"""
    transform_U(config::Matrix, index1::Int, index2::Int, sim::Simulation)

Builds diatomic Cartesian to internal coordinate transformation matrix as described in the SI of `10.1021/jacsau.0c00066`
"""
function transform_U(config::Matrix, index1::Int, index2::Int, sim::Simulation)
    masses = Structure.fractional_mass(sim, index1, index2)
    config[:, index2] .+=
        Structure.minimum_distance_translation(config, index1, index2, sim) # PBC wrap positions for correct transformation.
    r = transform_r(config, index1, index2)
    r1 = transform_r1(config, index1, index2)
    unity = LinearAlgebra.I(3)
    U_i1 = vcat(
        mapslices(x -> (x[1] - x[2]) / r, config[:, [index1, index2]]; dims = 2),
        mapslices(x -> (x[2] - x[1]) / r, config[:, [index2, index1]]; dims = 2),
    )
    U_i2 = vcat(
        mapslices(
            x -> (x[2] - x[1]) * (config[3, index2] - config[3, index1]) / (r^2 * r1),
            config[1:2, [index1, index2]];
            dims = 2,
        ),
        -r1 / r^2,
        mapslices(
            x -> (x[1] - x[2]) * (config[3, index2] - config[3, index1]) / (r^2 * r1),
            config[1:2, [index2, index1]];
            dims = 2,
        ),
        r1 / r^2,
    )
    U_i3 =
        [
            -(config[2, index1] - config[2, index2]),
            (config[1, index2] - config[1, index1]),
            0.0,
            -(config[2, index2] - config[2, index1]),
            (config[1, index1] - config[1, index2]),
            0.0,
        ] ./ (r1^2)
    return hcat(U_i1, U_i2, U_i3, vcat(unity .* masses[1], unity .* masses[2]))
end

export get_desorption_frame,
    get_desorption_angle,
    transform_from_internal_coordinates,
    transform_to_internal_coordinates,
    transform_U

end
