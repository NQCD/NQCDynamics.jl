export EANN_NOAu

"""
...
"""
struct EANN_NOAu <: FrictionModel
    path::String
    n_atoms::Int
    function EANN_NOAu(path::String, atoms::Atoms)
        initialize_NOAu_pes(path)
        new(path, length(atoms))
    end
end

function potential!(model::EANN_NOAu, V::AbstractVector, R::AbstractMatrix)
    coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
    cd(model.path) do
        calculate_energy_NOAu_pes!(coordinates_ang, reshape(V, (1,1)))
    end
end

function derivative!(model::EANN_NOAu, D::AbstractMatrix, R::AbstractMatrix)
    coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
    cd(model.path) do
        calculate_force_NOAu_pes!(coordinates_ang, zeros(1,1), D)
    end
    D .*= -1
end

function friction!(model::EANN_NOAu, F::AbstractMatrix, R::AbstractMatrix)
    #coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
    cd(model.path) do
        calculate_NOAu_friction_tensor!(R, F)
    end
end

function initialize_NOAu_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call init_pes function
        ccall((:init_pes_EFT_, "NOAu111_pes.dylib"), Cvoid, ())
    end
end

function deallocate_NOAu_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call deallocate_all function
        ccall((:deallocate_all_EFT_, "NOAu111_pes.dylib"), Cvoid, ())
    end
end

function calculate_energy_NOAu_pes!(coordinates::Array{Float64}, energy::Array{Float64})
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:get_energy_EFT_, "NOAu111_pes.dylib"),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}),
        coordinates, energy)
end

function calculate_force_NOAu_pes!(coordinates::Array{Float64}, energy::Array{Float64}, force::Array{Float64})
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:get_force_EFT_, "NOAu111_pes.dylib"),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
        coordinates, energy, force)
end

function calculate_NOAu_friction_tensor!(coordinates::AbstractMatrix, friction::AbstractMatrix)
    # EANN PES: calculate friction tensor
    ccall((:get_tensor_EFT_, "NOAu111_pes.dylib"),
        Cvoid,
        (Ref{Float64}, Ptr{Float64}),
        coordinates, friction)
end