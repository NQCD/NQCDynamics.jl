
module LDFA

using Dierckx
using DelimitedFiles
using DocStringExtensions
using ..Models
using PyCall

export CubeLDFAModel

const cube = PyNULL()

function __init__()
    pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)
    copy!(cube, pyimport("cube"))
end

abstract type LDFAModel <: Models.FrictionModel end

"""
$(TYPEDEF)

Wrapper for existing models that adds LDFA friction.
"""
struct CubeLDFAModel{T,M,S} <: LDFAModel
    model::M
    splines::S
    cube::PyObject
    ρ::Vector{T}
    radii::Vector{T}
    friction_atoms::Vector{Int}
end

function CubeLDFAModel(model::Models.Model, filename, atoms, friction_atoms=range(atoms))
    cube_object = cube.cube()
    cube_object.read(filename)

    ldfa_data, header = readdlm(joinpath(@__DIR__, "ldfa.txt"), ',', header=true)
    r = ldfa_data[:,1]
    splines = []
    for i in range(atoms)
        η = ldfa_data[:,atoms.numbers[i].+1]
        indices = η .!= ""
        ri = r[indices]
        η = η[indices]
        k = length(ri) > 3 ? 3 : length(ri) - 1
        push!(splines, Spline1D(ri, η; k=k))
    end

    ρ = zeros(length(atoms))
    radii = zero(ρ)

    CubeLDFAModel(model, splines, cube_object, ρ, radii, friction_atoms)
end

Models.potential!(model::LDFAModel, V::AbstractVector, R::AbstractMatrix) =
    Models.potential!(model.model, V, R)

Models.derivative!(model::LDFAModel, D::AbstractMatrix, R::AbstractMatrix) =
    Models.derivative!(model.model, D, R)

function density!(model::CubeLDFAModel, ρ::AbstractVector, R::AbstractMatrix)
    for (i, r) in enumerate(eachcol(R))
        ρ[i] = model.cube(r...)
    end
end

function Models.friction!(model::LDFAModel, F::AbstractMatrix, R::AbstractMatrix)
    density!(model, model.ρ, R)
    @. model.radii = 1 / (4/3 * π * model.ρ)^3
    DoFs = size(R, 1)
    for i in model.friction_atoms
        η = model.splines[i](model.radii[i])
        for j in axes(R, 1)
            F[(i-1)*DoFs+j, (i-1)*DoFs+j] = η
        end
    end
end

end # module