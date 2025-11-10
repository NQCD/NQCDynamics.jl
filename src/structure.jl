module Structure

using NQCDynamics: AbstractSimulation, Simulation, get_positions, get_velocities
using NQCBase
using RecursiveArrayTools: ArrayPartition
using ComponentArrays: ComponentVector
using LinearAlgebra
using Unitful
using UnitfulAtomic

int_or_index=Union{Int, CartesianIndex}

"""
    distance(config::Matrix, i1, i2)

Interatomic distance in Angstrom for a position Matrix.
"""
function distance(config::Matrix, i1::int_or_index, i2::int_or_index)
    return @views au_to_ang(norm(config[:,i1].-config[:,i2]))*u"Å"
end

"""
    distance(config::AbstractVector, i1::int_or_index, i2::int_or_index)
    Interatomic distance in Angstrom for DynamicsVariables. 
"""
function distance(config::AbstractVector, i1::int_or_index, i2::int_or_index)
    @views pos=NQCDynamics.get_positions(config)
    return au_to_ang(norm(pos[:,i1].-pos[:,i2]))*u"Å"
end

function distance(v1::Matrix, v2::Matrix)
    return au_to_ang(norm(v1-v2))*u"Å"
end

"""
    fractional_mass(sim::Simulation, index1::Int, index2::Int)

Returns m1/(m1+m2) and m2/(m1+m2) as a vector. 
"""
function fractional_mass(atoms::Atoms, index1::int_or_index, index2::int_or_index)
    return @views atoms.masses[[index1, index2]]./sum(atoms.masses[[index1, index2]])
end


"""
    reduced_mass(sim::Simulation, index1::Int, index2::Int)

Returns the reduced mass of the diatomic in atomic unit. 
"""
function reduced_mass(atoms::Atoms, index1::Int, index2::Int)
    return @views (sim.atoms.masses[index1]*sim.atoms.masses[index2])/sum(sim.atoms.masses[[index1, index2]])
end

"""
    minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, simulation::AbstractSimulation)

Outputs a translation vector to move config[:,ind2] such that the closest distance between `ind1` and `ind2` is reached. 
**The search of neighbouring unit cells will expand until the cutoff. If configurations are already subject to a minimum image convention, `cutoff=1` reduces unnecessary overhead. **
"""
function minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell;cutoff::Int=50)
    @debug "Cutoff was set to $(cutoff)"
    if cutoff==0 || isa(cell, InfiniteCell)
        return [0.0,0.0,0.0]
    end
    function generate_translations(order::Int)
        return Iterators.map(x->[x...],Iterators.product(-order:order, -order:order, -order:order))
    end
    min_distance=[]
    distance_converged=false
    cells=1
    while !(distance_converged || cells==cutoff+1)
        translations=generate_translations(cells)
        distances=map(x->norm(view(config, :,ind2)-view(config, :, ind1)+cell.Z*x), translations)
        push!(min_distance, (min(distances...), argmin(distances)))
        @views distance_converged=length(min_distance)>1 ? min_distance[end][1]≈min_distance[end-1][1] : false
        cells+=1
    end
    return cell.Z*collect(generate_translations(cells-1))[min_distance[end][2]]
end


"""
    angle_between(v1::Vector, v2::Vector)

Returns the angle between two vectors **in º** based on the scalar product definition. 
"""
function angle_between(v1::Vector, v2::Vector)
    return rad2deg(acos(dot(v1,v2)/prod(norm, [v1, v2])))
end


"""
    pbc_distance(config::Matrix, ind1::int_or_index, ind2::int_or_index, simulation::NQCDynamics.AbstractSimulation)

**Returns in Angstom, not in Bohr - Check units.**

Calculates the distance between two atoms, including a check if the copy of the second atom in any neighbouring unit cell is closer. 
This variant is designed for trajectories where cell boundary wrapping has been used. 
"""
function pbc_distance(config::Matrix, ind1::int_or_index, ind2::int_or_index, cell::PeriodicCell;args...)
    return @views auconvert(u"Å",norm(config[:, ind1]-config[:,ind2]-minimum_distance_translation(config, ind1, ind2, cell;args...)))
end

"""
    pbc_center_of_mass(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell, atoms::Atoms;args...)

Generates center of mass coordinates for two atoms, including a check if the copy of the second atom in any neighbouring unit cell is closer. 
"""
function pbc_center_of_mass(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell, atoms::Atoms;args...)
    pbc_translation=minimum_distance_translation(config, ind1, ind2, cell;args...)
    return @views (atoms.masses[ind1]*config[:,ind1]+atoms.masses[ind2]*(config[:,ind2]+pbc_translation))/sum(atoms.masses[[ind1,ind2]])
end

"""
    center_of_mass(config::Matrix, ind1::Int, ind2::Int, atoms::Atoms)

Generates center of mass coordinates for two atoms. 
"""
function center_of_mass(config::Matrix, ind1::Int, ind2::Int, atoms::Atoms)
    return @views (atoms.masses[ind1]*config[:,ind1]+atoms.masses[ind2]*config[:,ind2])/sum(atoms.masses[[ind1,ind2]])
end

"""
    velocity_center_of_mass(config::Matrix, ind1::Int, ind2::Int, simulation::NQCDynamics.AbstractSimulation)

`sum(m_i*v_i)/sum(m_i)`
"""
function velocity_center_of_mass(config::Matrix, ind1::Int, ind2::Int, atoms::Atoms)
    return center_of_mass(config, ind1, ind2, atoms)
end



"""
    fractional_mass(sim::NQCDynamics.AbstractSimulation, index1::Int, index2::Int)

Returns m1/(m1+m2) and m2/(m1+m2) as a vector. 
"""
function fractional_mass(sim::AbstractSimulation, index1::Int, index2::Int)
    return fractional_mass(sim.atoms,index1,index2)
end

function reduced_mass(sim::AbstractSimulation, index1::Int, index2::Int)
    return reduced_mass(sim.atoms,index1,index2)
end

function minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, simulation::AbstractSimulation;cutoff::Int=50)
    return minimum_distance_translation(config,ind1,ind2,simulation.cell;cutoff=cutoff)
end

function minimum_distance_translation(config::AbstractVector, ind1::Int, ind2::Int, simulation::AbstractSimulation;cutoff::Int=50)
    return minimum_distance_translation(get_positions(config),ind1,ind2,simulation.cell;cutoff=cutoff)
end

function pbc_distance(config::Matrix, ind1::int_or_index, ind2::int_or_index, sim::AbstractSimulation; args...)
    return pbc_distance(config,ind1,ind2,sim.cell; args...)
end

function pbc_distance(config::AbstractVector, ind1::int_or_index, ind2::int_or_index, sim::AbstractSimulation; args...)
    return pbc_distance(get_positions(config),ind1,ind2,sim.cell; args...)
end

function pbc_center_of_mass(config::Matrix, ind1::int_or_index, ind2::int_or_index, sim::AbstractSimulation; args...)
    return pbc_center_of_mass(config,ind1,ind2,sim.cell, sim.atoms; args...)
end

function pbc_center_of_mass(config::AbstractVector, ind1::int_or_index, ind2::int_or_index, sim::AbstractSimulation; args...)
    return pbc_center_of_mass(get_positions(config),ind1,ind2,sim.cell, sim.atoms; args...)
end

function velocity_center_of_mass(config::Matrix, ind1::int_or_index, ind2::int_or_index, sim::AbstractSimulation)
    return velocity_center_of_mass(config,ind1,ind2, sim.atoms)
end

function velocity_center_of_mass(config::AbstractVector, ind1::int_or_index, ind2::int_or_index, sim::AbstractSimulation)
    return velocity_center_of_mass(get_velocities(config),ind1,ind2, sim.atoms)
end


export minimum_distance_translation, pbc_distance, distance, angle_between, pbc_center_of_mass, velocity_center_of_mass, center_of_mass, reduced_mass, fractional_mass


end
