module NQCBaseCUDAext
using NQCBase
using CUDA
using Unitful, UnitfulAtomic
using PeriodicTable

# NQCBase.Atoms methods with CUDADevice

function NQCBase.Atoms(atom_types::AbstractVector{Symbol}, ::Type{NQCBase.CUDADevice})
    types = Vector{Symbol}(atom_types)
    numbers = Vector{Int}([element.number for element in elements[atom_types]])
    masses = Vector{Float64}([austrip(element.atomic_mass) for element in elements[atom_types]])
    return Atoms(types |> cu, numbers |> cu, masses |> cu)
end

function NQCBase.Atoms(masses::AbstractVector, ::Type{NQCBase.CUDADevice})
    natoms = length(masses)
    types = Vector{Symbol}(fill(:X, natoms))
    numbers = Vector{Int}(zeros(natoms))
    masses = Vector{Float64}(austrip.(masses))
    Atoms(types |> cu, numbers |> cu, masses |> cu)
end

end
