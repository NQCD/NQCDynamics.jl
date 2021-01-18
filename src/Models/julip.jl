using JuLIP
using PeriodicTable

export JuLIPModel

struct JuLIPModel{T} <: AdiabaticModel
    atoms::JuLIP.Atoms{T}
    function JuLIPModel(atoms::NonadiabaticMolecularDynamics.Atoms{N,T}, cell::PeriodicCell,
                        calculator::JuLIP.AbstractCalculator) where {N,T}
        jatoms = JuLIP.Atoms{T}(;
            X=zeros(3,length(atoms)),
            M=ustrip(auconvert.(u"u", atoms.masses)),
            Z=JuLIP.AtomicNumber.([elements[t].number for t in atoms.types]),
            cell=cell.vectors',
            pbc=cell.periodicity,
            calc=calculator)
        new{T}(jatoms)
    end
end

function potential!(model::JuLIPModel, V::AbstractVector, R::AbstractMatrix)
    set_positions!(model.atoms, R)
    V .= energy(model.atoms)
end

function derivative!(model::JuLIPModel, D::AbstractMatrix, R::AbstractMatrix)
    set_positions!(model.atoms, R)
    f = forces(model.atoms)
    for i=1:length(f)
        D[:,i] .= -f[i]
    end
end