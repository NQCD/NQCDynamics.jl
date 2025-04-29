"""

Add-on functions for NQCDistributions to make generating distributions more convenient. 

"""

import Distributions

"""
    VelocityBoltzmann(temperature, sim::AbstractSimulation; center = zeros(size(sim)))


Generates a VelocityBoltzmann distribution covering an entire Simulation.
If `NQCModels.mobileatoms` is modified, velocities for immobile atoms will always be zero.
"""
function NQCDistributions.VelocityBoltzmann(temperature, sim::AbstractSimulation; center=zeros(Float64, size(sim)))
    if length(NQCModels.mobileatoms(sim, 1)) == size(sim)[2]
        # Simple case: All atoms should be moving, so give them all Boltzmann distributions.
        return NQCDistributions.VelocityBoltzmann(temperature, sim.atoms.masses, size(sim); center=center)
    else
        # Not all atoms should be moving, so assign Boltzmann distribution only to mobile atoms.
        sampleable_array = convert(Matrix{Distributions.UnivariateDistribution}, hcat([[Distributions.Dirac(center[i, j]) for i in 1:size(sim)[1]] for j in 1:size(sim)[2]]...))
        for i in NQCModels.mobileatoms(sim, 1)
            for j in NQCModels.dofs(sim)
                sampleable_array[j, i] = NQCDistributions.VelocityBoltzmann(temperature, sim.atoms.masses[i]; center=center[j, i])
            end
        end
        return NQCDistributions.UnivariateArray(sampleable_array)
    end
end

function NQCDistributions.DynamicalDistribution(velocities, positions, simulation::Simulation; classical=Int[])
    mobile_atoms = NQCModels.mobileatoms(simulation.calculator.model, 1)
    simulation_dims = size(simulation)
    frozen_indices = symdiff(1:simulation_dims[2], mobile_atoms)
    if isempty(classical)
        return NQCDistributions.DynamicalDistribution(velocities, positions, simulation_dims; frozen_atoms=frozen_indices)
    else
        return NQCDistributions.DynamicalDistribution(velocities, positions, simulation_dims; frozen_atoms=frozen_indices, classical = classical)
    end

end
