using NonadiabaticMolecularDynamics
using PyCall
using DelimitedFiles
using Unitful
using UnitfulAtomic
using Plots
using DelimitedFiles
using Interpolations

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

data = readdlm(model_path * "/Tel_ag.txt")
time = 0:0.0025:9.9975
temperature = data[:,2]
fit = CubicSplineInterpolation(time, temperature)

function temperature_function(t)
    t_ps = ustrip(auconvert(u"ps", t))
    austrip(fit(t_ps)u"K")
end

println("Initialize...")
# model = EANN_H₂Ag(model_path, atoms)
model = FrictionHarmonic()
set_periodicity!(cell, [true, true, false])
sim = Simulation{TwoTemperatureMDEF}(atoms, model, temperature_function; cell=cell)
z = ClassicalDynamicals(vel_0, positions)
tspan = (0.0, 0.02u"ps")

@time solution = Dynamics.run_trajectory(z, tspan, sim; dt=0.05u"fs")

write_trajectory("trajectory.xyz", cell, atoms, get_positions.(solution.u))

ke = evaluate_kinetic_energy.(Ref(sim.atoms.masses), get_velocities.(solution.u))
pe = evaluate_potential_energy.(Ref(sim), get_positions.(solution.u))
plot(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", ke)), label="kinetic")
plot!(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", pe)), label="potential")
plot!(ustrip.(auconvert.(u"ps", solution.t)), ustrip.(auconvert.(u"eV", ke.+pe)), label="total")
