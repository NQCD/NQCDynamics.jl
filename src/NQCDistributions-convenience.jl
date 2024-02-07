"""
Add-on functions for NQCDistributions to make generating distributions more convenient. 

"""


"""
    VelocityBoltzmann(temperature, sim::AbstractSimulation; center = zeros(size(sim)))

Generates a VelocityBoltzmann distribution covering an entire Simulation. 
If `NQCModels.mobileatoms` is modified, velocities for immobile atoms will always be zero. 
"""
function VelocityBoltzmann(temperature, sim::AbstractSimulation; center = zeros(size(sim)))
	if mobileatoms(sim)!=size(sim,2)
		# Simple case: All atoms should be moving, so give them all Boltzmann distributions. 
		return VelocityBoltzmann(temperature, sim.atoms.masses, size(sim); center=center)
	else
		# Not all atoms should be moving, so assign Boltzmann distribution only to mobile atoms. 
		sampleable_array=convert(Matrix{UnivariateDistribution}, hcat([[Dirac(0.0) for i in 1:size(sim)[1]] for j in 1:size(sim)[2]]...))
		for i in mobileatoms(sim)
			sampleable_array[:,i].=VelocityBoltzmann(temperature, sim.atoms.masses[i]; center=center[:,i])
		end
		return sampleable_array
	end
end

export VelocityBoltzmann