using NonadiabaticMolecularDynamics
using Plots
using PeriodicTable
using UnitfulAtomic
using Unitful
using Random
using DelimitedFiles
using Distributions
# Random.seed!(10)
N = 150
mass = austrip(elements[:H].atomic_mass)
D = austrip(2u"eV")
a = austrip(1u"Å")
B = sqrt(D / 4π)*50
β = D
model = Models.Scattering1D(N=N, D=D, B=B, β=β, α=2β, a=a)

tspan = (0.0, 1000.0)

ke = austrip(10u"eV")
k = sqrt(2mass*ke)
# k = 30
position = Normal(3, 0)

function eval_loss(ensemble)
    avg = 0.0
    for e in ensemble
        # @assert e.retcode == :Terminated
        kinetic = 0.5 .* mass.* [get_velocities(u)[1] for u in e.u].^2
        ratio = kinetic[end] / kinetic[1]
        avg += ratio
    end
    avg / length(ensemble)
end

boundary(u,t,integrator) = (get_positions(u)[1] > 3) && (t > 100)
terminator = Dynamics.create_terminating_callback(boundary)

ks = 10:1:65
energies = []
for k in ks
    velocity = Normal(-k/mass, 0)
    distribution = InitialConditions.DynamicalDistribution(velocity, position, (1, 1))
    sim = Simulation{MDEF}(Atoms([:H]), model; DoFs=1)
    z = ClassicalDynamicals(rand(distribution)...)
    @time ensemble = run_ensemble(distribution, tspan, sim; dt=1, trajectories=1, callback=terminator)
    loss = eval_loss(ensemble)
    push!(energies, loss)
end

# open("friction_B=50B_dt1.txt", "w") do io
#     writedlm(io, [ks energies])
# end

plot(ustrip.(auconvert.(u"eV", ks.^2 ./2mass)), energies)
