"""
    QuantisedUtilities

    Holds common routines for QuantisedAtomic and QuantisedDiatomic

Inspired by VENUS96: [Hase1996](@cite)
"""
module QuantisedUtilities

"""
    separate_slab_and_molecule(atom_indices, r)

Get the coordinates of the molecule and slab separately.
"""
function separate_slab_and_molecule(atom_indices, r)
    molecule = r[:,atom_indices]
    slab_indices = [i for i in axes(r, 2) if i ∉ atom_indices]
    slab = r[:,slab_indices]
    return (molecule, slab)
end

function position_above_surface!(r, height, cell::PeriodicCell)
    r[3,:] .+= height
    a1 = cell.vectors[:,1]
    a2 = cell.vectors[:,2]
    displacement = rand()*a1+rand()*a2
    r .+= displacement
    return nothing
end


function position_above_surface!(r, height, ::InfiniteCell)
    r[3,:] .+= height
    return nothing
end

function apply_translational_impulse!(v, masses, translational_energy, direction)
    velocity_impulse = velocity_from_energy(masses, translational_energy)
    v .+= velocity_impulse .* normalize!(direction)
end

function velocity_from_energy(masses, energy)
    m = sum(masses)
    sqrt(2austrip(energy) / m)
end


end # module