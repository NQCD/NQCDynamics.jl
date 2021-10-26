using Test
using NonadiabaticMolecularDynamics
using ComponentArrays

atoms = Atoms([:C])
model = NonadiabaticModels.Harmonic()
cell = PeriodicCell(hcat(1))
sim = Simulation(atoms, model; cell=cell)

z = ComponentVector(v=fill(1.0, 1, length(atoms)), r=zeros(1, length(atoms)))

solution = run_trajectory(z, (0.0, 300), sim; dt=0.1)
Rs = vcat(get_positions.(solution.u)...)
@test !isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # check atom leaves cell

solution = run_trajectory(z, (0.0, 300), sim; callback=DynamicsUtils.CellBoundaryCallback(), dt=0.1)
Rs = vcat(get_positions.(solution.u)...)
@test isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # Check atom remains inside
