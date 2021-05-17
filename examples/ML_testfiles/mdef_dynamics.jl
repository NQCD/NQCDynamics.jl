using NonadiabaticMolecularDynamics
using Unitful
using PyCall

model_path = "examples/ML_testfiles"
friction_path = "examples/ML_testfiles/friction/"

cell, atoms, positions = read_system("examples/ML_testfiles/friction/H2Ag.xyz")
model = Models.SchNetPackModels.FrictionSchNetPackModel(model_path, friction_path, cell, atoms, [1,2])

sim = Simulation{MDEF}(atoms, model; cell=cell, temperature=500u"K")
z = ClassicalDynamicals(zero(positions), positions)

tspan = (0.0, 10.0) .* u"fs"
@time solution = Dynamics.run_trajectory(z, tspan, sim, dt=1u"fs")

write_trajectory("trajectory.xyz", cell, atoms, get_positions.(solution.u))
