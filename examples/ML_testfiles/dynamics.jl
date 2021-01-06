using NonadiabaticMolecularDynamics.IO
using NonadiabaticMolecularDynamics

model_path = "examples/ML_testfiles"

cell, atoms, positions = read_system("examples/ML_testfiles/H2.xyz")
model = Models.SchNetPackModel(model_path, cell, atoms)

sim = Simulation(atoms, model, Dynamics.Classical(); cell=cell)
z = Phasespace(positions, zero(positions))

tspan = (0.0, 1000.0)
@time solution = Dynamics.run_trajectory(z, tspan, sim)

write_trajectory("trajectory.xyz", cell, atoms, get_positions.(solution.u))
