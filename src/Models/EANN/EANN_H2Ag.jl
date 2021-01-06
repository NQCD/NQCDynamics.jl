export EANN_H₂Ag

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""

struct EANN_H₂Ag <: FrictionModel

    n_states::UInt
    potential!::Function
    derivative!::Function
    friction!::Function

    function EANN_H₂Ag(path::String, atoms::Atoms)

        initialize_H2Ag_pes(path)
        n_atoms = convert(Int, length(atoms))

        function potential!(V::AbstractVector, R::AbstractMatrix)
            coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
            cd(path) do
                calculate_H2Ag_pes_potential!(n_atoms, coordinates_ang, reshape(V, (1,1)))
            end
        end
        
        function derivative!(D::AbstractMatrix, R::AbstractMatrix)
            coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
            cd(path) do
                calculate_H2Ag_pes_forces!(n_atoms, coordinates_ang, D)
            end
            D .*= -1
        end

        function friction!(F::AbstractMatrix, R::AbstractMatrix)
            #coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
            cd(path) do
                calculate_H2Ag_friction_tensor!(R, F)
            end
        end
        
        new(1, potential!, derivative!, friction!)
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
