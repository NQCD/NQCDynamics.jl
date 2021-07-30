using Test
using NonadiabaticMolecularDynamics
using StatsBase: mean

ats = [3, 1, 10]
Ds = [1, 3, 7]
Ts = [0.5, 1, 10]

for (natoms, DoFs, T) in zip(ats, Ds, Ts)
	sim = Simulation(Atoms(rand(natoms)), Harmonic(); DoFs=DoFs, temperature=T)
	R0 = rand(DoFs, natoms)
	chain = InitialConditions.sample_configurations(sim, R0, 1e6)
	potential = evaluate_potential_energy.(sim, chain)
	@test mean(potential) / (DoFs*natoms) ≈ T / 2 rtol=1e-1
end
