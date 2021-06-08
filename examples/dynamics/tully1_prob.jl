using NonadiabaticMolecularDynamics
using Plots
using Unitful
using UnitfulAtomic
using Distributions
using LinearAlgebra

model = Models.TullyModelOne2()
atoms = Atoms([:H])
r = Normal(-10)
output_scattering_fssh = zeros(2, 2, 51)
output_scattering_ehrenfest = zeros(2, 2, 51)
ix = 0
mass = 2000

for k in 5:0.1:10
    print(k)
    global ix += 1
    sim_fssh = Simulation{FSSH}(atoms, model; DoFs=1)
    output_fssh = Ensembles.OutputScatteringFSSH()
    v = k/mass#sim_fssh.atoms.masses[1]
    distribution_fssh = InitialConditions.DynamicalDistribution(v, r, (1,1); state=1, type=:diabatic)
    selection_fssh = Ensembles.RandomSelection(distribution_fssh)
    reduction_fssh = Ensembles.MeanReduction()
    solution_fssh = Ensembles.run_ensemble(sim_fssh, (0.0, 3000.0), selection_fssh; trajectories=1e3,
        output=output_fssh, reduction=reduction_fssh, saveat=10.0,dtmax=5)
    print(" ")
    output_scattering_fssh[:,:,ix] = solution_fssh
 
    sim_ehrenfest = Simulation{Ehrenfest}(atoms, model; DoFs=1)
    output_ehrenfest = Ensembles.OutputScatteringEhrenfest()
    distribution_ehrenfest = InitialConditions.DynamicalDistribution(v, r, (1,1); state=1, type=:diabatic)
    selection_ehrenfest = Ensembles.RandomSelection(distribution_ehrenfest)
    reduction_ehrenfest = Ensembles.MeanReduction()
    solution_ehrenfest = Ensembles.run_ensemble(sim_ehrenfest, (0.0, 3000.0), selection_ehrenfest; trajectories=1e3,
        output=output_ehrenfest, reduction=reduction_ehrenfest, saveat=10.0,dtmax=5)
    output_scattering_ehrenfest[:,:,ix] = solution_ehrenfest
end

rng=5:0.1:10

plot_transmition1 = plot(rng, output_scattering_ehrenfest[1,2,:], title="Transmission 1", label="Ehrenfest", legend=:topleft)
plot!(rng, output_scattering_fssh[1,2,:], label="FSSH")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_scattering1 = plot(rng, output_scattering_ehrenfest[1,1,:], title="Scattering 1", label="Ehrenfest", legend=:bottomleft)
plot!(rng, output_scattering_fssh[1,1,:], label="FSSH")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_transmition2 = plot(rng, output_scattering_ehrenfest[2,2,:], title="Transmission 2", label="Ehrenfest", legend=:topleft)
plot!(rng, output_scattering_fssh[2,2,:], label="FSSH")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_scattering2 = plot(rng, output_scattering_ehrenfest[2,1,:], title="Scattering 2", label="Ehrenfest", legend=:topleft)
plot!(rng, output_scattering_fssh[2,1,:], label="FSSH")
ylabel!("Probability")
xlabel!("k / a.u.")

plot(plot_transmition1, plot_scattering1, plot_transmition2, plot_scattering2, layout = (2, 2), size=(1200,600))
