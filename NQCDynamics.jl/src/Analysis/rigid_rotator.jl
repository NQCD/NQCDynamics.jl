module RigidRotator
using Unitful, UnitfulAtomic
using NQCDynamics: Structure, AbstractSimulation
using LinearAlgebra: norm

"""
    classical_translational_energy(config::Any, ind1::Union{Int, CartesianIndex}, ind2::Union{Int, CartesianIndex}, sim::Simulation)

Returns the classical translational energy in Hartree
"""
function classical_translational_energy(config::Any, ind1::Union{Int,CartesianIndex}, ind2::Union{Int,CartesianIndex}, sim :: AbstractSimulation)
    return sum(sim.atoms.masses[[ind1, ind2]]) / 2 * norm(velocity_center_of_mass(config, ind1, ind2, sim.atoms))^2
end

"""
    classical_rotation_energy(J::Union{Int, CartesianIndex}, config::Any, ind1::Union{Int, CartesianIndex}, ind2::Union{Int, CartesianIndex}, sim::Simulation)

Classical rotation energy of a rigid diatomic rotor in Hartree
"""
function classical_rotation_energy(J::Int, config::Any, ind1::Union{Int,CartesianIndex}, ind2::Union{Int,CartesianIndex}, sim :: AbstractSimulation)
    I = reduced_mass(sim, ind1, ind2) * austrip(pbc_distance(config, ind1, ind2, sim))^2 # Classical moment of inertia in atomic units.
    return J * (J + 1) / (2 * I)
end

"""
    harmonic_vibration_energy(ν::Union{Int, CartesianIndex}, k::Float, ind1::Union{Int, CartesianIndex}=1, ind2::Union{Int, CartesianIndex}=2, sim::Simulation)

Vibrational energy of a harmonic oscillator with the force constant k and vibrational level ν.
"""
function harmonic_vibration_energy(ν::Union{Int}, k::Float64, ind1::Union{Int,CartesianIndex}, ind2::Union{Int,CartesianIndex}, sim :: AbstractSimulation)
    return (ν + 1 / 2) * sqrt(k / reduced_mass(sim, ind1, ind2))
end

export classical_translational_energy, classical_rotation_energy, harmonic_vibration_energy

end
