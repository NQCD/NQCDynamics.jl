using NonadiabaticMolecularDynamics
using PyCall
using DelimitedFiles
using Unitful
using UnitfulAtomic
using Plots

vel_0 = zeros(18,3)
ps = 98.22694788464062
model_path = "/Users/wojciechstark/Desktop/h2ag_dyn_files"
input = model_path * "/start.xyz"
velo_path = model_path * "/init.200"
cell, atoms, positions = read_system(input)

#read initial velocities
data = open(readdlm, velo_path)
vel_0[1,:] = data[1, 5:7]
vel_0[2,:] = data[2, 5:7]
vel_0 = austrip.(vel_0.*u"Å/ps")
vel_0 = copy(transpose(vel_0))
vel_0[1,:] .*= atoms.masses[1]
vel_0[2,:] .*= atoms.masses[2]

println("Initialize...")
cb, vals = Dynamics.create_energy_saving_callback()
model = Models.EANN_H₂Ag(model_path, atoms)
sim = Simulation(atoms, model, Dynamics.MDEF{Float64}(length(atoms)); temperature=0u"K", cell=cell)
z = Phasespace(positions, vel_0)
tspan = (0.0, 5.0u"ps")

@time solution = Dynamics.run_trajectory(z, tspan, sim; callback=cb, dt=0.05u"fs", reltol=1e-8)

write_trajectory("trajectory.xyz", cell, atoms, get_positions.(solution.u))
plot(solution.t, evaluate_hamiltonian.(Ref(sim), solution.u))
ke = NonadiabaticMolecularDynamics.evaluate_kinetic_energy.(Ref(sim), get_momenta.(solution.u))
pe = NonadiabaticMolecularDynamics.evaluate_configurational_energy.(Ref(sim), get_positions.(solution.u))
plot(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", ke)), label="kinetic")
plot!(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", pe)), label="potential")
plot!(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", ke.+pe)), label="total")
