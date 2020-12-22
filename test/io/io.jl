using Test
using NonadiabaticMolecularDynamics.IO

cell, atoms, R = read_system("io/slab.xyz")
@test size(R)[1] == 3
@test atoms isa Atoms
@test cell isa PeriodicCell
