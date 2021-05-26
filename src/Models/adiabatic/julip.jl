import .JuLIP
using PeriodicTable
using StaticArrays

export JuLIPModel

"""
    JuLIPModel(atoms::Atoms{N,T}, cell::PeriodicCell,
               calculator::JuLIP.AbstractCalculator) where {N,T}

Interface to [`JuLIP.jl`](https://github.com/JuliaMolSim/JuLIP.jl).
This might be broken due to updates in JuLIP.

```julia
import JuLIP
at = JuLIP.bulk(:Si, cubic=true)

R = ang_to_au.(Matrix(hcat(at.X...)))
atoms = Atoms(JuLIP.chemical_symbol.(at.Z))
cell = PeriodicCell(ang_to_au.(at.cell))
calculator = JuLIP.StillingerWeber()

model = Models.JuLIPModel(atoms, cell, calculator)
```
"""
mutable struct JuLIPModel{T,A<:NamedTuple,B<:NamedTuple} <: AdiabaticModel
    atoms::JuLIP.Atoms{T}
    tmp::A
    tmp_d::B
    function JuLIPModel(atoms::NonadiabaticMolecularDynamics.Atoms{N,T}, cell::PeriodicCell,
                        calculator::JuLIP.AbstractCalculator) where {N,T}

        jatoms = JuLIP.Atoms{T}(;
            X=zeros(3,length(atoms)),
            P=zeros(3,length(atoms)),
            M=au_to_u.(atoms.masses),
            Z=JuLIP.AtomicNumber.([elements[t].number for t in atoms.types]),
            cell=au_to_ang.(cell.vectors'),
            pbc=cell.periodicity,
            calc=calculator)

        tmp = JuLIP.alloc_temp(calculator, jatoms)
        tmp_d = JuLIP.alloc_temp_d(calculator, jatoms)

        new{T,typeof(tmp),typeof(tmp_d)}(jatoms, tmp, tmp_d)
    end
end

function potential!(model::JuLIPModel, V::AbstractVector, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    try
        V[1] = JuLIP.energy!(model.tmp, model.atoms.calc, model.atoms)
    catch e
        if e isa BoundsError
            model.tmp = JuLIP.alloc_temp(model.atoms.calc, model.atoms)
            V[1] = JuLIP.energy!(model.tmp, model.atoms.calc, model.atoms)
        end
    end
    V[1] = eV_to_au(V[1])
    return V
end

function derivative!(model::JuLIPModel, D::AbstractMatrix, R::AbstractMatrix)
    JuLIP.set_positions!(model.atoms, au_to_ang.(R))
    try
        JuLIP.forces!(JuLIP.vecs(D), model.tmp_d, model.atoms.calc, model.atoms)
    catch e
        if e isa BoundsError
            model.tmp_d = JuLIP.alloc_temp_d(model.atoms.calc, model.atoms)
            JuLIP.forces!(JuLIP.vecs(D), model.tmp_d, model.atoms.calc, model.atoms)
        end
    end
    D .= -eV_per_ang_to_au.(D)
    return D
end
