
using NQCBase: au_to_ang, eV_to_au, eV_per_ang_to_au
using .JuLIP: JuLIP

export JuLIPModel

"""
    struct JuLIPModel{T} <: ClassicalModel

Model for interfacing with JuLIP potentials.
"""
struct JuLIPModel{T} <: ClassicalModel
    atoms::JuLIP.Atoms{T}
end

NQCModels.ndofs(::JuLIPModel) = 3

function NQCModels.potential(model::JuLIPModel, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    return eV_to_au(JuLIP.energy(model.atoms))
end

function NQCModels.potential!(model::JuLIPModel, V::Matrix{<:Number}, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    V .= eV_to_au(JuLIP.energy(model.atoms))
end

function NQCModels.derivative!(model::JuLIPModel, D::AbstractMatrix, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    force = JuLIP.forces(model.atoms)
    D .= -eV_per_ang_to_au.(JuLIP.mat(force))
    return D
end
