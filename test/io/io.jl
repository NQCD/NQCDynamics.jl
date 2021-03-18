using Test
using NonadiabaticMolecularDynamics
using PyCall
using TypedTables

path = joinpath(dirname(pathof(NonadiabaticMolecularDynamics)), "../test")
cell, atoms, R = read_system(joinpath(path, "io", "slab.xyz"))
@test size(R)[1] == 3
@test atoms isa Atoms
@test cell isa PeriodicCell

tab = Table((velocity=[1, 23, 4, 5, 6], position=[1, 4, 5, 56, 67]))
out = SystemStore(tab, atoms, cell)

variables = InitialConditions.DynamicalDistribution(5, 10, (3,2))
out = SystemStore(variables, atoms, cell)

save("test.jld2", out)
loaded = load("test.jld2")
@test out.data == loaded.data
@test loaded.atoms == atoms
rm("test.jld2")