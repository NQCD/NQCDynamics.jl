using NonadiabaticMolecularDynamics

atoms = Atoms{Float64}([:H, :C, :Pt])
model = Models.Harmonic()
sim = Simulation(atoms, model, Dynamics.Classical())

test_motion!(sim)