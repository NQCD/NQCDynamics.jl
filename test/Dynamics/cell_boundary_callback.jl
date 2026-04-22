using Test
using NQCDynamics
using ComponentArrays
import DiffEqBase

atoms = Atoms([:C])
model = NQCModels.Harmonic()
cell = PeriodicCell(hcat(1))
sim = Simulation(atoms, model; cell=cell)

z = ComponentVector(v=fill(1.0, 1, length(atoms)), r=zeros(1, length(atoms)))

solution = run_dynamics(sim, (0.0, 300), z; output=OutputPosition, dt=0.1)
Rs = reduce(vcat, solution[:OutputPosition])
@test !isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # check atom leaves cell

solution = run_dynamics(sim, (0.0, 300), z; output=OutputPosition, callback=DiffEqBase.CallbackSet(DynamicsUtils.CellBoundaryCallback(),), dt=0.1)
Rs = reduce(vcat, solution[:OutputPosition])
@assert length(Rs) == length(solution[:OutputPosition]) # Check callback doesn't mess with output length
@test isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # Check atom remains inside
