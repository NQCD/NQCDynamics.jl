export EANN_H₂Ag

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""
struct EANN_H₂Ag{S,T} <: FrictionModel
    path::String
    atoms::Atoms{S,T}
    h2indices::Vector{UInt}
    tmp_coordinates::Matrix{T}
    tmp_friction_coordinates::Matrix{T}
    tmp_friction::Matrix{T}
    tmp_energy::Matrix{T}
    function EANN_H₂Ag(path::String, atoms::Atoms{S,T}) where {S,T}
        initialize_H2Ag_pes(path)
        h2indices = findall(atoms.types .== :H)
        @assert h2indices == [1, 2]
        new{S,T}(path, atoms, h2indices,
                 zeros(3, length(atoms)),
                 zeros(3, 2),
                 zeros(6, 6),
                 zeros(1, 1))
    end
end

function potential!(model::EANN_H₂Ag, V::AbstractVector, R::AbstractMatrix)
    model.tmp_coordinates .= au_to_ang.(R)
    cd(model.path) do
        calculate_H2Ag_pes_potential!(length(model.atoms), model.tmp_coordinates, model.tmp_energy)
    end
    V[1] = eV_to_au(model.tmp_energy[1])
end

function derivative!(model::EANN_H₂Ag, D::AbstractMatrix, R::AbstractMatrix)
    model.tmp_coordinates .= au_to_ang.(R)
    cd(model.path) do
        calculate_H2Ag_pes_forces!(length(model.atoms), model.tmp_coordinates, D)
    end
    D .= eV_per_ang_to_au.(D)
end

function friction!(model::EANN_H₂Ag, F::AbstractMatrix, R::AbstractMatrix)
    @views model.tmp_friction_coordinates .= au_to_ang.(R[:,model.h2indices])
    cd(model.path) do
        calculate_H2Ag_friction_tensor!(model.tmp_friction_coordinates, model.tmp_friction)
    end
    F[1:6,1:6] .= model.atoms.masses[1] .* model.tmp_friction
    F .= ps_inv_to_au.(F)
end

function initialize_H2Ag_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call init_pes function
        ccall((:pes_init_, "h2ag111_pes.dylib"), Cvoid, ())
    end
end

function calculate_H2Ag_pes_potential!(n_atoms::Int64, coordinates::Matrix{Float64}, energy::Matrix{Float64})
    # EANN PES: get energy from EANN to vars: energy and force
    ccall((:pot0_, "h2ag111_pes.dylib"),
        Cvoid,
        (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
        n_atoms, coordinates, energy)
end

function calculate_H2Ag_pes_forces!(n_atoms::Int64, coordinates::Matrix{Float64}, force::Matrix{Float64})
    # EANN PES: get force from EANN to vars: energy and force
    ccall((:dpeshon_, "h2ag111_pes.dylib"),
        Cvoid,
        (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
        n_atoms, coordinates, force)
end

function calculate_H2Ag_friction_tensor!(coordinates::AbstractMatrix, friction::AbstractMatrix)
    # EANN PES: calculate friction tensor
    ccall((:tensor_, "h2ag111_pes_friction.dylib"),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}),
        coordinates, friction)
end
