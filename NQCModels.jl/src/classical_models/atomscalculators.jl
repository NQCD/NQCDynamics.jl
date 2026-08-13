using AtomsCalculators
import AtomsBase
using NQCBase
using Unitful, UnitfulAtomic
using StaticArrays

# NQCD uses these units
AtomsCalculators.energy_unit(M::ClassicalModel,) = u"hartree"
AtomsCalculators.length_unit(M::ClassicalModel,) = u"a0_au"

struct AtomsCalculatorsModel{C} <: ClassicalModel
    calc_object::C
    energy_unit
    length_unit
    ndofs::Int
    atoms::NQCBase.Atoms
    cell::NQCBase.AbstractCell
end

"""
    AtomsCalculatorsModel(calc_object, structure::AtomsBase.AbstractSystem)

Model interface to AtomsCalculators. Supply this function with the calculator object you would use with AtomsCalculators and the correct unit conversions will be automatically applied.
"""
function AtomsCalculatorsModel(calc_object, structure::AtomsBase.AbstractSystem)
    atoms = NQCBase.Atoms(structure)
    cell = NQCBase.Cell(structure)
    return AtomsCalculatorsModel(
        calc_object,
        AtomsCalculators.energy_unit(calc_object),
        AtomsCalculators.length_unit(calc_object),
        3, # AtomsBase is always 3D
        atoms,
        cell,
    )
end

function AtomsCalculatorsModel(calc_object, structure::NQCBase.Structure)
    atoms = structure.atoms
    cell = structure.cell
    return AtomsCalculatorsModel(
        calc_object,
        AtomsCalculators.energy_unit(calc_object),
        AtomsCalculators.length_unit(calc_object),
        3, # AtomsBase is always 3D
        atoms,
        cell,
    )
end

NQCModels.ndofs(::AtomsCalculatorsModel) = 3

function NQCModels.potential(model::AtomsCalculatorsModel, R::AbstractMatrix)
    # Convert into system format expected by AtomsBase
    sy = NQCBase.System(model.atoms, R, model.cell)
    # Convert AtomsBase calculator energy unit back out.
    return austrip(AtomsCalculators.potential_energy(sy, model.calc_object))
end

function NQCModels.potential!(model::AtomsCalculatorsModel, V::Matrix{<:Number}, R::AbstractMatrix)
    V .= NQCModels.potential(model, R)
end

function NQCModels.derivative!(model::AtomsCalculatorsModel, D::AbstractMatrix, R::AbstractMatrix)
    # Convert into system format expected by AtomsBase
    sy = NQCBase.System(model.atoms, R, model.cell)
    forces = .-reduce(hcat, AtomsCalculators.forces(sy, model.calc_object)) # Convert to matrix representation rather than Vector{Vector}.
    D .= austrip.(forces) # Convert back into atomic units.
    return D
end


# Minimal AtomsBase Calculator implementation for Classical Models (so we don't have to worry about which state to select. )
function AtomsCalculators.potential_energy(
    sys::AtomsBase.AbstractSystem,
    model::ClassicalModel,
)
    nqcd_pos = NQCBase.Position(sys)
    return NQCModels.potential(model, nqcd_pos) * AtomsCalculators.energy_unit(model)
end

function AtomsCalculators.forces(
    sys::AtomsBase.AbstractSystem,
    model::ClassicalModel,
)
    nqcd_pos = NQCBase.Position(sys)
    forces_with_unit = -NQCModels.derivative(model, nqcd_pos) .* AtomsCalculators.energy_unit(model) ./ AtomsCalculators.length_unit(model)
    return SVector{3}.(eachcol(forces_with_unit))
end

function AtomsCalculators.virial(
    sys::AtomsBase.AbstractSystem,
    model::ClassicalModel,
)
    nd::Int = NQCModels.ndofs(model)
    return zeros(SMatrix{nd,nd})
end
