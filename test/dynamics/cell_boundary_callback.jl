using Test
using NonadiabaticMolecularDynamics

atoms = Atoms([:C])
model = Models.Harmonic()
cell = PeriodicCell(hcat(1))
sim = Simulation(atoms, model, Classical(); DoFs=1, cell=cell)

z = ClassicalDynamicals(fill(1.0, 1, length(atoms)), zeros(1, length(atoms)))

solution = Dynamics.run_trajectory(z, (0.0, 300), sim; dt=0.1)
Rs = vcat(get_positions.(solution.u)...)
@test !isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # check atom leaves cell

solution = Dynamics.run_trajectory(z, (0.0, 300), sim; callback=CellBoundaryCallback, dt=0.1)
Rs = vcat(get_positions.(solution.u)...)
@test isapprox(minimum(Rs), 0, atol=1e-4) | (minimum(Rs) >= 0) # Check atom remains inside
