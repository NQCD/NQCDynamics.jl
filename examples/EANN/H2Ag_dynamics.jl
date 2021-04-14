using NonadiabaticMolecularDynamics
using PyCall
using DelimitedFiles
using Unitful
using UnitfulAtomic
using Plots

vel_0 = zeros(18,3)
model_path = "."
input = model_path * "/start.xyz"
velo_path = model_path * "/init.200"
cell, atoms, positions = read_system(input)

#read initial velocities
data = open(readdlm, velo_path)
vel_0[1,:] = data[1, 5:7]
vel_0[2,:] = data[2, 5:7]
vel_0 = austrip.(vel_0.*u"Å/ps")
vel_0 = copy(transpose(vel_0))

println("Initialize...")
model = Models.EANN_H₂Ag(model_path, atoms)
sim = Simulation{MDEF}(atoms, model, cell=cell)
z = ClassicalDynamicals(vel_0, positions)
tspan = (0.0, 0.1u"ps")

@time solution = Dynamics.run_trajectory(z, tspan, sim; dt=5u"fs")

write_trajectory("trajectory.xyz", cell, atoms, get_positions.(solution.u))

ke = evaluate_kinetic_energy.(Ref(sim.atoms.masses), get_velocities.(solution.u))
pe = evaluate_potential_energy.(Ref(sim), get_positions.(solution.u))
plot(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", ke)), label="kinetic")
plot!(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", pe)), label="potential")
plot!(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", ke.+pe)), label="total")
