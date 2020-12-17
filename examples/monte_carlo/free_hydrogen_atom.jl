#=
    Normal mode sampling of a free ring polymer
    
Here we have a single hydrogen atom fixed in place and
perform PIMC to sample the normal modes.
=#
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using NonadiabaticMolecularDynamics.IO
using Unitful
using Plots

model = Models.Analytic.Free()
atoms = Atoms.AtomicParameters(Atoms.InfiniteCell(), [:H])

Δ = Dict([(:H, 1.0)]) # This is for the centroid mode and does nothing here as the centroid is fixed.
beads = 100
sys = RingPolymerSystem{MonteCarlo}(atoms, model, 30u"K", beads, Δ; passes=1000, fix=[1])

positions = reshape(zeros(3), 3, 1) # Create zero coordinates
R = cat([positions for i=1:beads]..., dims=3) # Copy coordinates for every bead
output = InitialConditions.run_monte_carlo_sampling(sys, R)

@show output.acceptance
write_trajectory("sampling.xyz", output.R, atoms)
plot(output.energy)