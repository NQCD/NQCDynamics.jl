#=
    Normal mode sampling of a free ring polymer
    
Here we have a single hydrogen atom fixed in place and
perform PIMC to sample the normal modes.
=#
using NonadiabaticMolecularDynamics
using Unitful
using Plots
using PyCall

model = Models.Free()
atoms = Atoms{Float64}([:H])

Δ = Dict([(:H, 1.0)]) # This is for the centroid mode and does nothing here as the centroid is fixed.
beads = 10
monte_carlo = InitialConditions.PathIntegralMonteCarlo{Float64}(Δ, length(atoms), 300, [1], 0.5, beads)
sim = RingPolymerSimulation(atoms, model, Dynamics.Classical(), beads; temperature=30u"K")

positions = reshape(zeros(3), 3, 1) # Create zero coordinates
R = cat([positions for i=1:beads]..., dims=3) # Copy coordinates for every bead
output = InitialConditions.run_monte_carlo_sampling(sim, monte_carlo, R)

@show output.acceptance
write_trajectory("sampling.xyz", InfiniteCell(), atoms, output.R)
plot(output.energy)