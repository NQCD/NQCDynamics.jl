#push!(LOAD_PATH, pwd())

using NonadiabaticMolecularDynamics.IO
#using NonadiabaticMolecularDynamics.Dynamics
using NonadiabaticMolecularDynamics.Systems
using NonadiabaticMolecularDynamics.Models.EANN_H2_Cu
#using DifferentialEquations
using Unitful
using UnitfulAtomic


model_path = "/Users/wojciechstark/Desktop/H2_on_Cu/1_h2cu_pes"
input_f = model_path * "/H2Cu_example.xyz"
atoms, positions = read_system(input_f)
model = EANN_H2_Cu.EannH2CuModel(model_path, atoms)
p = Systems.System(atoms, model)
positions_ang = copy(ustrip(auconvert.(u"Å", positions)))

print("\nPositions: \n")
print(positions_ang)
print("\nEnergy: \n")
en = model.get_V0(positions_ang)
print(en)
print("\nForces: \n")
f = model.get_D0(positions_ang)
print(f)


#tspan = (0.0, 2000.0)
#problem = ODEProblem(Dynamics.differential!, phasespace, tspan, p)
#solution = solve(problem)
#write_trajectory("trajectory.xyz", solution, atoms)
