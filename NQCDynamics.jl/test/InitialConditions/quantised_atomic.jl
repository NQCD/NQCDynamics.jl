using Test
using NQCDynamics
using NQCDynamics.InitialConditions.ConfigureAtomic
using LinearAlgebra: norm
using Unitful
using UnitfulAtomic


atoms = Atoms([:O])
model = Free()
sim = Simulation(atoms, model)
surface = ConfigureAtomic.SurfaceParameters(sim.atoms.masses, [1], Matrix{Float64}(undef, 3, 0), 10.0, [0, 0, 1.0])