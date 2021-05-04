export EANN_H₂Ag
using PeriodicTable
using Libdl

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""
struct EANN_H₂Ag{S,T} <: AdiabaticFrictionModel
    potential_function::Ptr{Nothing}
    force_function::Ptr{Nothing}
    friction_function::Ptr{Nothing}
    potential_path::String
    friction_path::String
    atoms::Atoms{S,T}
    h2indices::Vector{UInt}
    tmp_coordinates::Matrix{T}
    tmp_friction_coordinates::Matrix{T}
    tmp_friction::Matrix{T}
    tmp_energy::Matrix{T}
    function EANN_H₂Ag(potential_path::String, friction_path::String, atoms::Atoms{S,T}) where {S,T}
        pes_lib = dlopen(potential_path)
        friction_lib = dlopen(friction_path)

        initialize!(potential_path, pes_lib)

        potential_function = dlsym(pes_lib, "pot0_")
        force_function = dlsym(pes_lib, "dpeshon_")
        friction_function = dlsym(friction_lib, "tensor_")
    
        h2indices = findall(atoms.types .== :H)
        @assert h2indices == [1, 2]
        new{S,T}(potential_function, force_function, friction_function,
                 potential_path, friction_path, atoms, h2indices,
                 zeros(3, 2),
                 zeros(3, 2),
                 zeros(6, 6),
                 zeros(1, 1))
    end
end

function initialize!(path, pes_lib)
    cd(splitdir(path)[1]) do
        # EANN PES: call init_pes function
        pes_init = dlsym(pes_lib, "pes_init_")
        ccall(pes_init, Cvoid, ())
    end
end

function potential!(model::EANN_H₂Ag, V::AbstractVector, R::AbstractMatrix)
    @views model.tmp_coordinates .= au_to_ang.(R[:,model.h2indices])
    cd(splitdir(model.potential_path)[1]) do
        ccall(model.potential_function, Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
              2, model.tmp_coordinates, model.tmp_energy)
    end
    V[1] = eV_to_au(model.tmp_energy[1])
end

function derivative!(model::EANN_H₂Ag, D::AbstractMatrix, R::AbstractMatrix)
    @views model.tmp_coordinates .= au_to_ang.(R[:,model.h2indices])
    cd(splitdir(model.potential_path)[1]) do
        ccall(model.force_function, Cvoid, (Ref{Int64}, Ref{Float64}, Ptr{Float64}),
              2, model.tmp_coordinates, D)
    end
    D .= eV_per_ang_to_au.(D)
end

function friction!(model::EANN_H₂Ag, F::AbstractMatrix, R::AbstractMatrix)
    @views model.tmp_friction_coordinates .= au_to_ang.(R[:,model.h2indices])
    cd(splitdir(model.friction_path)[1]) do
        ccall(model.friction_function, Cvoid, (Ref{Float64}, Ptr{Float64}),
              model.tmp_coordinates, model.tmp_friction)
    end
    F[1:6,1:6] .= ps_inv_to_au.(model.tmp_friction)
    F .*= austrip(elements[:H].atomic_mass)
end

