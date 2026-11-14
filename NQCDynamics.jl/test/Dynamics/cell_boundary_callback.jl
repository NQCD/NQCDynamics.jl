using Test
using NQCDynamics
using ComponentArrays

atoms = Atoms([:C])
model = NQCModels.Harmonic()
cell = PeriodicCell(hcat(1))
sim = Simulation(atoms, model; cell=cell)

z = ComponentVector(v=fill(1.0, 1, length(atoms)), r=zeros(1, length(atoms)))

solution = run_dynamics(sim, (0.0, 300), z; output=OutputPosition, dt=0.1)
Rs = reduce(vcat, solution[:OutputPosition])
@test !isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # check atom leaves cell

solution = run_dynamics(sim, (0.0, 300), z; output=OutputPosition, callback=DynamicsUtils.CellBoundaryCallback(), dt=0.1)
Rs = reduce(vcat, solution[:OutputPosition])
@test isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # Check atom remains inside
