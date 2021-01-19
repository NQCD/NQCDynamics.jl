using Test
using NonadiabaticMolecularDynamics
using PyCall

path = joinpath(dirname(pathof(NonadiabaticMolecularDynamics)), "../test")
cell, atoms, R = read_system(joinpath(path, "io", "slab.xyz"))
@test size(R)[1] == 3
@test atoms isa Atoms
@test cell isa PeriodicCell
