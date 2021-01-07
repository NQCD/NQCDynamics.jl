export EANN_H₂Ag

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""
struct EANN_H₂Ag <: FrictionModel
    path::String
    n_atoms::Int
    function EANN_H₂Ag(path::String, atoms::Atoms)
        initialize_H2Ag_pes(path)
        new(path, length(atoms))
    end
end

function potential!(model::EANN_H₂Ag, V::AbstractVector, R::AbstractMatrix)
    coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
    cd(model.path) do
        calculate_H2Ag_pes_potential!(model.n_atoms, coordinates_ang, reshape(V, (1,1)))
    end
end

function derivative!(model::EANN_H₂Ag, D::AbstractMatrix, R::AbstractMatrix)
    coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
    cd(model.path) do
        calculate_H2Ag_pes_forces!(model.n_atoms, coordinates_ang, D)
    end
    D .*= -1
end

function friction!(model::EANN_H₂Ag, F::AbstractMatrix, R::AbstractMatrix)
    #coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
    cd(model.path) do
        calculate_H2Ag_friction_tensor!(R, F)
    end
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
