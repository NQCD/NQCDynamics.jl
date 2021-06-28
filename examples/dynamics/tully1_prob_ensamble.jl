using NonadiabaticMolecularDynamics
using Plots
using Unitful
using UnitfulAtomic
using Distributions
using LinearAlgebra
using ProgressMeter

model = TullyModelOne()
atoms = Atoms(2000)
r = Normal(-4)
reduction = Ensembles.MeanReduction()
n_beads=5
n_traj = 2000
momenta = 3:0.5:30

output_scattering_fssh = zeros(2, 2, length(momenta))
output_scattering_ehrenfest = zeros(2, 2, length(momenta))
output_scattering_erpmd = zeros(2, 2, length(momenta))

@showprogress 1 "Running ensembles..." for (i, k) in enumerate(momenta)
    e1=(k^2)/(2*atoms.masses[1])
    v = k / atoms.masses[1]
    distribution = InitialConditions.DynamicalDistribution(v, r, (1,1); state=1, type=:diabatic)
    distribution2 = InitialConditions.DynamicalDistribution(v, r, (1,1,n_beads); state=1, type=:diabatic)
    selection = Ensembles.RandomSelection(distribution)
    selection2 = Ensembles.RandomSelection(distribution2)
    distance = 10
    time = distance / v
    tspan = (0.0, time)

    sim_fssh = Simulation{FSSH}(atoms, model; DoFs=1)
    output = Ensembles.OutputStateResolvedScattering1D(sim_fssh)
    solution_fssh = Ensembles.run_ensemble(sim_fssh, tspan, selection;
        trajectories=n_traj, output=output, reduction=reduction)
    output_scattering_fssh[:,:,i] .= solution_fssh.u
 
    sim_ehrenfest = Simulation{Ehrenfest}(atoms, model; DoFs=1)
    output = Ensembles.OutputStateResolvedScattering1D(sim_ehrenfest)
    solution_ehrenfest = Ensembles.run_ensemble(sim_ehrenfest, tspan, selection;
        trajectories=n_traj, output=output, reduction=reduction)
    output_scattering_ehrenfest[:,:,i] .= solution_ehrenfest.u

    sim_erpmd = RingPolymerSimulation{Ehrenfest}(atoms, model, n_beads; DoFs=1, temperature=10u"K")
    output = Ensembles.OutputStateResolvedScattering1D(sim_erpmd)
    solution_erpmd = Ensembles.run_ensemble(sim_erpmd, tspan, selection2;
        trajectories=n_traj, output=output, reduction=reduction)
    output_scattering_erpmd[:,:,i] .= solution_erpmd.u
end

plot_transmition1 = plot(momenta, output_scattering_ehrenfest[1,2,:], title="Transmission 1", label="Ehrenfest", legend=:topright)
plot!(momenta, output_scattering_fssh[1,2,:], label="FSSH")
plot!(momenta, output_scattering_erpmd[1,2,:], label="Ehrenfest-RPMD")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_scattering1 = plot(momenta, output_scattering_ehrenfest[1,1,:], title="Scattering 1", label="Ehrenfest", legend=:topright)
plot!(momenta, output_scattering_fssh[1,1,:], label="FSSH")
plot!(momenta, output_scattering_erpmd[1,1,:], label="Ehrenfest-RPMD")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_transmition2 = plot(momenta, output_scattering_ehrenfest[2,2,:], title="Transmission 2", label="Ehrenfest", legend=:topleft)
plot!(momenta, output_scattering_fssh[2,2,:], label="FSSH")
plot!(momenta, output_scattering_erpmd[2,2,:], label="Ehrenfest-RPMD")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_scattering2 = plot(momenta, output_scattering_ehrenfest[2,1,:], title="Scattering 2", label="Ehrenfest", legend=:topright)
plot!(momenta, output_scattering_fssh[2,1,:], label="FSSH")
plot!(momenta, output_scattering_erpmd[2,1,:], label="Ehrenfest-RPMD")
ylabel!("Probability")
xlabel!("k / a.u.")

plot(plot_transmition1, plot_scattering1, plot_transmition2, plot_scattering2, layout = (2, 2), size=(1200,600))
