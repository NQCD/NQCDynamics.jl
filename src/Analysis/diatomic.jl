"""
Analysis functions for surface chemistry of diatomic molecules. 
"""
module Diatomic
using NQCDynamics: AbstractSimulation, Simulation, get_positions, Structure
using NQCBase
using Unitful, UnitfulAtomic
using LinearAlgebra

"""
    get_desorption_frame(trajectory::Vector, diatomic_indices::Vector{Int}, simulation::AbstractSimulation;surface_normal::Vector=[0,0,1], surface_distance_threshold=5.0*u"Å")

Determines the index in a trajectory where surface desorption begins. 

This is evaluated using two conditions:

1. In the trajectory, the diatomic must be `surface_distance_threshold` or further away from the highest other atom. (In `surface_normal` direction). 

2. Desorption begins at the turning point of the centre of mass velocity component along `surface_normal`, indicating overall movement away from the surface. 
"""
function get_desorption_frame(trajectory::Vector, diatomic_indices::Vector{Int}, simulation::AbstractSimulation;surface_normal::Vector = [0,0,1], surface_distance_threshold = 5.0*u"Å")
    # Two criteria: Distance between two atoms (wrt PBC) must be below `distance_threshold` and CoM velocity must be positive wrt surface normal vector. 
    function surface_distance_condition(x)
        highest_z = max(get_positions(x)[3, symdiff(1:end, diatomic_indices)]...)
        if abs(au_to_ang(Structure.pbc_center_of_mass(x, diatomic_indices..., simulation)[3] - highest_z)) * u"Å" ≥ surface_distance_threshold
            @debug "Surface distance condition evaluated true"
            return true
        else
            return false
        end
    end
    function com_velocity_condition(x)
        if dot(Structure.velocity_center_of_mass(x, diatomic_indices..., simulation), normalize(surface_normal)) > 0
            @debug "Normal velocity condition evaluated true"
            return true
        else
            return false
        end
    end
    desorbed_frame = findfirst(surface_distance_condition, trajectory)
    if isa(desorbed_frame, Nothing)
        @debug "No desorption event found. "
        return nothing
    else
        @debug "Found desorption event in frame $(desorbed_frame)"
        leaving_surface_frame = findfirst(!com_velocity_condition, reverse(trajectory[1:desorbed_frame-1]))
        if isa(leaving_surface_frame, Nothing)
            @debug "Couldn't find a sub-zero CoM velocity component relative to surface normal."
            return nothing
        else
            return desorbed_frame - leaving_surface_frame
        end
    end
end

function get_desorption_angle(trajectory::AbstractVector, indices::Vector{Int}, simulation::AbstractSimulation; surface_normal = [0,0,1],  surface_distance_threshold = 5.0*u"Å")
    # First determine where the reaction occurred on the surface. 
    desorption_frame = get_desorption_frame(trajectory, indices, simulation;surface_distance_threshold = surface_distance_threshold, surface_normal = surface_normal)
    if isa(desorption_frame, Nothing)
        @debug "No desorption event detected in trajectory"
        return nothing
    end
    @debug "Desorption frame: $(desorption_frame)"
    # Determine the average centre of mass velocity to decrease error due to vibration and rotation orthogonal to true translational component. 
    com_velocities = zeros(eltype(trajectory[1]), length(surface_normal), length(trajectory)-desorption_frame)
    for i in 1:length(trajectory)-desorption_frame
        com_velocities[:, i] .= velocity_center_of_mass(trajectory[i+desorption_frame], indices[1], indices[2], simulation)
    end
    average_velocity=mean(com_velocities;dims=2)
    # Now convert into an angle by arccos((a•b)/(|a|*|b|))
    return Structure.angle_between(vec(average_velocity), surface_normal)
end

# 6✖6 Cartesian to internal coordinate transformation. 
function transform_to_internal_coordinates(to_transform::Matrix, config::Matrix, index1::Int, index2::Int, sim::Simulation)
    U = transform_U(config, index1, index2, sim) # Generate transformation matrix
    return transpose(U) * to_transform * U
end

function transform_from_internal_coordinates(to_transform::Matrix, config::Matrix, index1::Int, index2::Int, sim::Simulation)
    U_inv = inv(transform_U(config, index1, index2, sim))
    return transpose(U_inv) * to_transform * U_inv
end

function transform_r(config, index1::Int, index2::Int)
    return sqrt(sum(mapslices(x->(x[2]-x[1])^2, config[:, [index1, index2]];dims=2)))
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
    config[:,index2] .+= Structure.minimum_distance_translation(config, index1, index2, sim) # PBC wrap positions for correct transformation. 
    r = transform_r(config, index1, index2)
    r1 = transform_r1(config, index1, index2)
    unity = LinearAlgebra.I(3)
    U_i1 = vcat(
        mapslices(x->(x[1]-x[2])/r, config[:, [index1, index2]]; dims=2),
        mapslices(x->(x[2]-x[1])/r, config[:, [index2, index1]]; dims=2),
    )
    U_i2 = vcat(
        mapslices(x->(x[2]-x[1])*(config[3, index2]-config[3, index1])/(r^2*r1), config[1:2, [index1, index2]]; dims=2),
        -r1/r^2,
        mapslices(x->(x[1]-x[2])*(config[3, index2]-config[3, index1])/(r^2*r1), config[1:2, [index2, index1]]; dims=2),
        r1/r^2,
    )
    U_i3 = [
        -(config[2, index1]-config[2, index2]),
        (config[1, index2]-config[1, index1]),
        0.0,
        -(config[2, index2]-config[2, index1]),
        (config[1, index1]-config[1, index2]),
        0.0,
    ] ./ (r1^2)
    return hcat(U_i1, U_i2, U_i3, vcat(unity.*masses[1], unity.*masses[2]))
end

export get_desorption_frame, get_desorption_angle, transform_from_internal_coordinates, transform_to_internal_coordinates, transform_U

end