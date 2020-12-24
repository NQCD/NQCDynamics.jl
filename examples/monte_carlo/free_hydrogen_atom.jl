#=
    Normal mode sampling of a free ring polymer
    
Here we have a single hydrogen atom fixed in place and
perform PIMC to sample the normal modes.
=#
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.IO
using Unitful
using Plots

model = Models.Free()
cell = InfiniteCell{Float64}()
atoms = Atoms{Float64}([:H])

Δ = Dict([(:H, 1.0)]) # This is for the centroid mode and does nothing here as the centroid is fixed.
beads = 100
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), 300, [1])
sim = RingPolymerSimulation(1, 30u"K", cell, atoms, model, monte_carlo, beads)

positions = reshape(zeros(3), 3, 1) # Create zero coordinates
R = cat([positions for i=1:beads]..., dims=3) # Copy coordinates for every bead
output = InitialConditions.run_monte_carlo_sampling(sim, R)

@show output.acceptance
write_trajectory("sampling.xyz", cell, atoms, output.R)
plot(output.energy)