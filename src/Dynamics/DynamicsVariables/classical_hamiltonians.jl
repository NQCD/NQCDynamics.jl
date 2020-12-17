export evaluate_configurational_energy
export get_spring_energy
export get_potential_energy

function evaluate_configurational_energy(system::System, R::Matrix)
    Electronics.evaluate_potential(system.model, R)
end

function evaluate_configurational_energy(system::RingPolymerSystem, R::Array{T, 3}) where {T}
    get_spring_energy(system, R) + get_potential_energy(system, R)
end

function get_potential_energy(system::RingPolymerSystem, R::Array{T, 3}) where {T}
    V = 0.0
    @views for i=1:n_beads(system)
        V += Electronics.evaluate_potential(system.model, R[:,:,i])
    end
    V
end

"""
    get_spring_energy(system::RingPolymerSystem, R::RingPolymerArray)
    
Calculate the ring polymer spring potential.

I don't like having to write these loops but it seems faster than any alternative.
"""
function get_spring_energy(system::RingPolymerSystem, R::Array{T, 3}) where {T}
    E = 0.0
    for bead=1:n_beads(system)-1
        for i in system.ring_polymer.quantum_atoms # Only for quantum nuclei
            for j=1:n_DoF(system)
                E += masses(system)[i] * (R[j, i, bead] - R[j, i, bead+1])^2
            end
        end
    end
    for i in system.ring_polymer.quantum_atoms
        for j=1:n_DoF(system)
            E += masses(system)[i] * (R[j, i, n_beads(system)] - R[j, i, 1])^2
        end
    end
    E * system.ring_polymer.ω_n^2/2
end
