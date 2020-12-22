export evaluate_configurational_energy
export get_spring_energy
export get_potential_energy

function evaluate_configurational_energy(sim::Simulation, R::Matrix)
    Calculators.evaluate_potential(sim.calculator, R)
end

function evaluate_configurational_energy(sim::RingPolymerSimulation, R::Array{T, 3}) where {T}
    get_spring_energy(sim, R) + get_potential_energy(sim, R)
end

function get_potential_energy(sim::RingPolymerSimulation, R::Array{T, 3}) where {T}
    V = 0.0
    @views for i=1:length(sim.beads)
        V += Calculators.evaluate_potential(sim.calculator, R[:,:,i])
    end
    V
end

"""
    get_spring_energy(system::RingPolymerSystem, R::RingPolymerArray)
    
Calculate the ring polymer spring potential.

I don't like having to write these loops but it seems faster than any alternative.
"""
function get_spring_energy(sim::RingPolymerSimulation, R::Array{T, 3}) where {T}
    E = 0.0
    for bead=1:length(sim.beads)-1
        for i in sim.beads.quantum_atoms # Only for quantum nuclei
            for j=1:sim.DoFs
                E += sim.atoms.masses[i] * (R[j, i, bead] - R[j, i, bead+1])^2
            end
        end
    end
    for i in sim.beads.quantum_atoms
        for j=1:sim.DoFs
            E += sim.atoms.masses[i] * (R[j, i, length(sim.beads)] - R[j, i, 1])^2
        end
    end
    E * sim.beads.ω_n^2/2
end
