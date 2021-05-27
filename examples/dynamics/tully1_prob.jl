using NonadiabaticMolecularDynamics
using Plots
using Unitful
using UnitfulAtomic
using Distributions
using LinearAlgebra

function scattering_output(sol, i)
    final= last(sol.u) # get final configuration from trajectory
    σ = get_density_matrix(final)
    populations = real.(diag(σ)) # Get adiabatic populations from diagonal of density matrix
    output = zeros(2, 2) # Initialise output, left column reflection on state 1/2, right column transmission
    x = get_positions(final)[1]
    if x > 0 # If final position past 0 then we count as transmission 
        output[:,2] .= populations
    else # If final position left of 0 then we count as reflection
        output[:,1] .= populations
    end
    return output #(output, false)
end

model = Models.TullyModelOne()
atoms = Atoms([:H])
r = Normal(-10)
output_scattering_fssh = zeros(2, 2, 301)
output_scattering_ehrenfest = zeros(2, 2, 301)
ix = 0

for k in 0:0.1:30
    print(k)
    global ix += 1
    sim_fssh = Simulation{FSSH}(atoms, model; DoFs=1)
    output_fssh = Ensembles.OutputDiabaticPopulation(sim_fssh)
    v = k/sim_fssh.atoms.masses[1]
    distribution_fssh = InitialConditions.DynamicalDistribution(v, r, (1,1); state=1, type=:diabatic)
    selection_fssh = Ensembles.RandomSelection(distribution_fssh)
    reduction_fssh = Ensembles.MeanReduction()
    solution_fssh = Ensembles.run_ensemble(sim_fssh, (0.0, 3000.0), selection_fssh; trajectories=2e3,
        reduction=reduction_fssh, saveat=10.0)
    output_scattering_fssh[:,:,ix] = scattering_output(solution_fssh, 0)

    sim_ehrenfest = Simulation{Ehrenfest}(atoms, model; DoFs=1)
    output_ehrenfest = Ensembles.OutputDiabaticPopulation(sim_ehrenfest)
    distribution_ehrenfest = InitialConditions.DynamicalDistribution(v, r, (1,1); state=1, type=:diabatic)
    selection_ehrenfest = Ensembles.RandomSelection(distribution_ehrenfest)
    reduction_ehrenfest = Ensembles.MeanReduction()
    solution_ehrenfest = Ensembles.run_ensemble(sim_ehrenfest, (0.0, 3000.0), selection_ehrenfest; trajectories=2e3,
        reduction=reduction_ehrenfest, saveat=10.0)
    output_scattering_ehrenfest[:,:,ix] = scattering_output(solution_ehrenfest, 0)
end

print(output_scattering_ehrenfest)
print(" HERE STARTS FSSH: ")
print(output_scattering_fssh)

plot_transmition = plot(0:0.1:30, output_scattering_ehrenfest[1,2,:], title="Transmission", label="Ehrenfest")
plot!(0:0.1:30, output_scattering_fssh[1,2,:], label="FSSH")
ylabel!("Probability")
xlabel!("k / a.u.")

plot_scattering = plot(0:0.1:30, output_scattering_ehrenfest[1,1,:], title="Scattering", label="Ehrenfest")
plot!(0:0.1:30, output_scattering_fssh[1,1,:], label="FSSH")
ylabel!("Probability")
xlabel!("k / a.u.")

plot(plot_transmition, plot_scattering, layout = (1, 2), size=(1200,400))
