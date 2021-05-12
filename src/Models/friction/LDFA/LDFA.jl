
"""
Models within attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""
module LDFA

using DataInterpolations
using DelimitedFiles
using DocStringExtensions
using UnitfulAtomic, Unitful
using ....NonadiabaticMolecularDynamics: PeriodicCell, au_to_ang, apply_cell_boundaries!
using ..Models
using ..PyCall

export CubeLDFAModel

const cube = PyNULL()

function __init__()
    pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)
    copy!(cube, pyimport("cube"))
end

abstract type LDFAModel <: Models.AdiabaticFrictionModel end

"""
$(TYPEDEF)

Wrapper for existing models that adds LDFA friction.

This model uses a cube file to evaluate the electron density used to calculate the friction.

$(FIELDS)
"""
struct CubeLDFAModel{T,M,S} <: LDFAModel
    """Model that provides the energies and forces."""
    model::M
    """Splines fitted to the numerical LDFA data."""
    splines::S
    """Cube object provided by `cube.py`."""
    cube::PyObject
    """Temporary array for storing the electron density."""
    ρ::Vector{T}
    """Temporary array for the Wigner-Seitz radii."""
    radii::Vector{T}
    """Indices of atoms that should have friction applied."""
    friction_atoms::Vector{Int}
    "Optional vector that specifies offset that should be applied."
    cube_offset::Vector{T}
    "Periodic cell containing the electron density."
    cell::PeriodicCell{T}
    function CubeLDFAModel(model, splines, cube, ρ, radii, friction_atoms, cube_offset, cell)
        new{eltype(cell),typeof(model),typeof(splines)}(
            model, splines, cube, ρ, radii, friction_atoms, cube_offset, cell)
    end
end

function CubeLDFAModel(model::Models.Model, filename, atoms, cell;
                       friction_atoms=collect(range(atoms)),
                       cube_offset=zeros(3).*u"Å")
    cube_object = cube.cube()
    cube_object.read(filename)

    ldfa_data, header = readdlm(joinpath(@__DIR__, "ldfa.txt"), ',', header=true)
    r = ldfa_data[:,1]
    splines = []
    for i in range(atoms)
        η = ldfa_data[:,atoms.numbers[i].+1]
        indices = η .!= ""
        ri = convert(Vector{Float64}, r[indices])
        η = convert(Vector{Float64}, η[indices])
        push!(ri, 10.0) # Ensure it goes through 0.0 for large r.
        push!(η, 0.0)
        push!(splines, CubicSpline(η, ri))
    end

    ρ = zeros(length(atoms))
    radii = zero(ρ)

    CubeLDFAModel(model, splines, cube_object, ρ, radii, friction_atoms, ustrip.(uconvert.(u"Å", cube_offset)), cell)
end

Models.potential!(model::LDFAModel, V::AbstractVector, R::AbstractMatrix) =
    Models.potential!(model.model, V, R)

Models.derivative!(model::LDFAModel, D::AbstractMatrix, R::AbstractMatrix) =
    Models.derivative!(model.model, D, R)

function density!(model::CubeLDFAModel, ρ::AbstractVector, R::AbstractMatrix)
    for i in model.friction_atoms
        r = R[:,i]
        apply_cell_boundaries!(model.cell, r)
        r .= au_to_ang.(r)
        ρ[i] = austrip(model.cube((r .+ model.cube_offset)...) * u"Å^-3")
    end
end

function Models.friction!(model::LDFAModel, F::AbstractMatrix, R::AbstractMatrix)
    density!(model, model.ρ, R)
    clamp!(model.ρ, 0, Inf)
    @. model.radii = 1 / cbrt(4/3 * π * model.ρ)
    DoFs = size(R, 1)
    for i in model.friction_atoms
        η = model.radii[i] < 10 ? model.splines[i](model.radii[i]) : 0.0
        for j in axes(R, 1)
            F[(i-1)*DoFs+j, (i-1)*DoFs+j] = η
        end
    end
    return F
end

end # module