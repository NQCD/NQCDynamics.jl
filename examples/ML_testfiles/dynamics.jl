using NonadiabaticMolecularDynamics.IO
using NonadiabaticMolecularDynamics.Dynamics
using NonadiabaticMolecularDynamics.Models.ML
using NonadiabaticMolecularDynamics.Systems
using DifferentialEquations

model_path = "examples/ML_testfiles"

atoms, phasespace = read_system("examples/ML_testfiles/H2.xyz")
model = SchNetPackModel(model_path, atoms)

p = Systems.System(atoms, model)

tspan = (0.0, 20000.0)
problem = ODEProblem(Dynamics.differential!, phasespace, tspan, p)
solution = solve(problem)

write_trajectory("trajectory.xyz", solution, atoms)