using Test
using NonadiabaticMolecularDynamics
using StaticArrays

atoms = Atoms{Float64}([:C, :H])
@test atoms isa Atoms{2, Float64}
@test length(atoms) == 2
@test atoms.numbers == [6, 1]
@test atoms.types == [:C, :H]
@test atoms.masses isa SVector

@test range(atoms) == 1:2
