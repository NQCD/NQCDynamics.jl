module Structure

using NQCDynamics: AbstractSimulation, Simulation, get_positions, get_velocities
using NQCBase
using LinearAlgebra
using Unitful
using UnitfulAtomic

int_or_index=Union{Int, CartesianIndex}

"""
    distance(config::Matrix, i1, i2)
    Interatomic distance in Angstrom for a position Matrix. 
TBW
"""
function distance(config::Matrix, i1::int_or_index, i2::int_or_index)
    return au_to_ang(norm(config[:,i1].-config[:,i2]))*u"Å"
end

function distance(v1::Matrix, v2::Matrix)
    return au_to_ang(norm(v1-v2))*u"Å"
end

"""
    fractional_mass(sim::Simulation, index1::Int, index2::Int)

Returns m1/(m1+m2) and m2/(m1+m2) as a vector. 
"""
function fractional_mass(atoms::Atoms, index1::int_or_index, index2::int_or_index)
    return atoms.masses[[index1, index2]]./sum(atoms.masses[[index1, index2]])
end


"""
    reduced_mass(sim::Simulation, index1::Int, index2::Int)

Returns the reduced mass of the diatomic in atomic unit. 
"""
function reduced_mass(atoms::Atoms, index1::Int, index2::Int)
    return (atoms.masses[index1]*atoms.masses[index2])/sum(atoms.masses[[index1, index2]])
end

#ToDo NQCDynamics aware version Simulation -> Atoms


#ToDo Variant with Simulation -> PeriodicCell
"""
    minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, simulation::NQCDynamics.AbstractSimulation)

Outputs a translation vector to move config[:,ind2] such that the closest distance between `ind1` and `ind2` is reached. 
**The search of neighbouring unit cells will expand until the cutoff. If configurations are already subject to a minimum image convention, `cutoff=1` reduces unnecessary overhead. **
"""
function minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell;cutoff::Int=50)
    @debug "Cutoff was set to $(cutoff)"
    if cutoff==0 || isa(cell, InfiniteCell)
        return [0.0,0.0,0.0]
    end
    translations=[]
    min_distance=[]
    distance_converged=false
    search=1
    while !(distance_converged || search==cutoff+1)
        translations=Iterators.map(x->[x...],Iterators.product(-search:search, -search:search, -search:search))
        distances=map(x->norm(config[:,ind2].-config[:,ind1]+cell.vectors*x), translations)
        push!(min_distance, (min(distances...), argmin(distances)))
        distance_converged=length(min_distance)>1 ? min_distance[end][1]≈min_distance[end-1][1] : false
        search+=1
    end
    @debug "PBC distance check ended after $(search-1) iterations."
    if search==cutoff+1
        @debug "PBC distance check terminated due to cutoff. If you aren't PBC-wrapping positions, check your cutoff is sufficiently high. "
    end
    @debug "Translation of choice is $(collect(translations)[min_distance[end][2]])"
    return cell.vectors*collect(translations)[min_distance[end][2]]
end

function minimum_distance_translation(config::Matrix, ind1::Int, ind2::Int, cell::InfiniteCell;cutoff::Int=50)
	config
end


"""
    angle_between(v1::Vector, v2::Vector)

Returns the angle between two vectors **in º** based on the scalar product definition. 
"""
function angle_between(v1::Vector, v2::Vector)
    return rad2deg(acos(dot(v1,v2)/prod(norm, [v1, v2])))
end

#ToDo Versions for Any --> Matrix
#ToDo Versions for Simulation --> Cell, Masses

"""
    pbc_distance(config::Matrix, ind1::Int, ind2::Int, simulation::NQCDynamics.AbstractSimulation)

**Returns in Angstom, not in Bohr - Check units.**

Calculates the distance between two atoms, including a check if the copy of the second atom in any neighbouring unit cell is closer. 
This variant is designed for trajectories where cell boundary wrapping has been used. 
"""
function pbc_distance(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell;args...)
    return auconvert(u"Å",norm(config[:, ind1]-config[:,ind2]-minimum_distance_translation(config, ind1, ind2, cell;args...)))
end

"""
    pbc_center_of_mass(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell, atoms::Atoms;args...)

Generates center of mass coordinates for two atoms, including a check if the copy of the second atom in any neighbouring unit cell is closer. 
"""
function pbc_center_of_mass(config::Matrix, ind1::Int, ind2::Int, cell::PeriodicCell, atoms::Atoms;args...)
    pbc_translation=minimum_distance_translation(config, ind1, ind2, cell;args...)
    return (atoms.masses[ind1]*config[:,ind1]+atoms.masses[ind2]*(config[:,ind2]+pbc_translation))/sum(atoms.masses[[ind1,ind2]])
end

"""
    center_of_mass(config::Matrix, ind1::Int, ind2::Int, atoms::Atoms)

Generates center of mass coordinates for two atoms. 
"""
function center_of_mass(config::Matrix, ind1::Int, ind2::Int, atoms::Atoms)
    return (atoms.masses[ind1]*config[:,ind1]+atoms.masses[ind2]*config[:,ind2])/sum(atoms.masses[[ind1,ind2]])
end

"""
    velocity_center_of_mass(config::Matrix, ind1::Int, ind2::Int, simulation::NQCDynamics.AbstractSimulation)

    `sum(m_i*v_i)/sum(m_i)`
TBW
"""
function velocity_center_of_mass(config::Matrix, ind1::Int, ind2::Int, atoms::Atoms)
    return center_of_mass(config, ind1, ind2, atoms)
end

struct OutputSubsetKineticEnergy
    indices::Vector{Int}
end
function (kinetic_energy::OutputSubsetKineticEnergy)(sol, i)
    return map(x->DynamicsUtils.classical_kinetic_energy(sol.prob.p.atoms.masses[kinetic_energy.indices], x), [DynamicsUtils.get_velocities(i) for i in sol.u])
end


export minimum_distance_translation, pbc_distance, distance, angle_between, pbc_center_of_mass, velocity_center_of_mass, center_of_mass, OutputSubsetKineticEnergy, reduced_mass, fractional_mass


end