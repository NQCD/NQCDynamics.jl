export evaluate_hamiltonian
export get_spring_energy
export evaluate_potential_energy
export evaluate_kinetic_energy

function evaluate_hamiltonian(sim::AbstractSimulation, u::DynamicalVariables)
    evaluate_hamiltonian(sim, get_velocities(u), get_positions(u))
end

function evaluate_hamiltonian(sim::AbstractSimulation, v::AbstractMatrix, r::AbstractMatrix)
    k = evaluate_kinetic_energy(sim.atoms.masses, v)
    e = evaluate_potential_energy(sim, r)
    k + e
end

function evaluate_hamiltonian(sim::RingPolymerSimulation, v::Array{T,3}, r::Array{T,3}) where {T}
    E = 0.0
    @views for i=1:length(sim.beads)
        E += evaluate_kinetic_energy(sim.atoms.masses, v[:,:,i])
    end
    E += evaluate_potential_energy(sim, r)
    E / length(sim.beads)
end

evaluate_kinetic_energy(masses, v) = sum(masses' .* v.^2)/2

function evaluate_potential_energy(sim::AbstractSimulation, R::AbstractMatrix)
    Calculators.evaluate_potential!(sim.calculator, R)
    sim.calculator.potential[1]
end

function evaluate_potential_energy(sim::RingPolymerSimulation, R::Array{T, 3}) where {T}
    Calculators.evaluate_potential!(sim.calculator, R)
    get_spring_energy(sim, R) + sum(sim.calculator.potential)[1]
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
