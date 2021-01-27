using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using DelimitedFiles
using Unitful
using Plots

vel_0 = zeros(18,3)
ps = 98.22694788464062
fs = 0.09822694788464063

model_path = "/Users/wojciechstark/Desktop/h2ag_dyn_files"
input = model_path * "/start.xyz"
velo_path = model_path * "/init.200"
cell, atoms, positions = read_system(input)

#read initial velocities
data = open(readdlm, velo_path)
vel_0[1,:] = data[1, 5:7]
vel_0[2,:] = data[2, 5:7]
vel_0 ./= ps
vel_0 = copy(transpose(vel_0))
vel_0[1,:] .*= atoms.masses
vel_0[2,:] .*= atoms.masses

println("Initialize...")
model = Models.EANN_H₂Ag(model_path, atoms)
sim = Simulation(atoms, model, Dynamics.MDEF{Float64}(length(atoms)); temperature=0u"K", cell=cell)
z = Phasespace(positions, vel_0)

tspan = (0.0, austrip(10u"fs"))
dt = austrip(0.05u"ps")

#@time solution, e_pot = Dynamics.run_trajectory(z, tspan, sim, dt)
@time solution = Dynamics.run_trajectory(z, tspan, sim)

write_trajectory("trajectory.xyz", cell, atoms, get_positions.(solution.u))
#plot(solution.t, e_pot)