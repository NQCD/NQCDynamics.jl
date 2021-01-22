import .JuLIP
using PeriodicTable

export JuLIPModel

struct JuLIPModel{T} <: AdiabaticModel
    atoms::JuLIP.Atoms{T}
    function JuLIPModel(atoms::NonadiabaticMolecularDynamics.Atoms{N,T}, cell::PeriodicCell,
                        calculator::JuLIP.AbstractCalculator) where {N,T}
        jatoms = JuLIP.Atoms{T}(;
            X=zeros(3,length(atoms)),
            M=ustrip.(auconvert.(u"u", atoms.masses)),
            Z=JuLIP.AtomicNumber.([elements[t].number for t in atoms.types]),
            cell=cell.vectors',
            pbc=cell.periodicity,
            calc=calculator)
        new{T}(jatoms)
    end
end

function potential!(model::JuLIPModel, V::AbstractVector, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, ustrip.(auconvert.(u"Å", R)))
    V .= austrip(JuLIP.energy(model.atoms)u"eV")
end

function derivative!(model::JuLIPModel, D::AbstractMatrix, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, ustrip.(auconvert.(u"Å", R)))
    f = JuLIP.forces(model.atoms)
    for i=1:length(f)
        D[:,i] .= -austrip.(f[i]u"eV/Å")
    end
end