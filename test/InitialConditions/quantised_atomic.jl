using Test
using NQCDynamics
using NQCDynamics.InitialConditions.QuantisedAtomic
using LinearAlgebra: norm
using Unitful
using UnitfulAtomic


atoms = Atoms([:O])
model = Free()
sim = Simulation(atoms, model)
surface = QuantisedAtomic.SurfaceParameters(sim.atoms.masses, [1], Matrix{Float64}(undef, 3, 0), 10.0, [0, 0, 1.0])