using Test
using NQCBase

atoms = Atoms([:C, :H])
@test atoms isa Atoms{Vector{Symbol}, Vector{Int}, Vector{Float64}}
@test length(atoms) == 2
@test atoms.numbers == [6, 1]
@test atoms.types == [:C, :H]
@test atoms.masses isa Vector

@test range(atoms) == 1:2

@test atoms[1] isa Atoms{Vector{Symbol}, Vector{Int}, Vector{Float64}}
@test atoms[1:2] isa Atoms{Vector{Symbol}, Vector{Int}, Vector{Float64}}

@test Atoms(2000) isa Atoms{Vector{Symbol}, Vector{Int}, Vector{Float64}}
