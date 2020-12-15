export evaluate_configurational_energy
export get_spring_energy
export get_potential_energy

function evaluate_configurational_energy(system::System, R::Matrix)
    Electronics.evaluate_potential(system.model, R)
end

function evaluate_configurational_energy(system::RingPolymerSystem, R::RingPolymerArray)
    get_spring_energy(system, R) + get_potential_energy(system, R)
end

function get_potential_energy(system::RingPolymerSystem, R::RingPolymerArray)
    V = 0.0
    for i=1:n_beads(system)
        V += Electronics.evaluate_potential(system.model, R.x[i])
    end
    V
end

"""
    get_spring_energy(system::RingPolymerSystem, R::RingPolymerArray)
    
Calculate the ring polymer spring potential.

I don't like having to write these loops but it seems faster than any alternative.
"""
function get_spring_energy(system::RingPolymerSystem, R::RingPolymerArray)
    E = 0.0
    for bead=1:n_beads(system)-1
        for i in system.ring_polymer.quantum_atoms # Only for quantum nuclei
            for j=1:n_DoF(system)
                E += masses(system)[i] * (R.x[bead][j, i] - R.x[bead+1][j, i])^2
            end
        end
    end
    for i in system.ring_polymer.quantum_atoms
        for j=1:n_DoF(system)
            E += masses(system)[i] * (R.x[n_beads(system)][j, i] - R.x[1][j, i])^2
        end
    end
    E * system.ring_polymer.ω_n^2/2
end
