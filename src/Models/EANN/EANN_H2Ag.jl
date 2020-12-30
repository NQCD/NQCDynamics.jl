export EannH2AgModel

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""

struct EannH2AgModel <: AdiabaticModel

    n_states::UInt
    potential!::Function
    derivative!::Function

    function EannH2AgModel(path::String, atoms::Atoms)

        n_atoms = convert(Int, length(atoms))
        
        initialize_H2Ag_pes(path)

        function potential!(V::Vector, R::AbstractMatrix)
            coordinates_ang = copy(ustrip(auconvert.(u"Å", R))) # converts coordinates into Angstrom
            cd(path) do
                calculate_H2Ag_pes_potential!(n_atoms, coordinates_ang, reshape(V, (1,1)))
            end
        end
        
        function derivative!(D::AbstractMatrix, R::AbstractMatrix)
            coordinates_ang = copy(ustrip(auconvert.(u"Å", R))) # converts coordinates into Angstrom
            cd(path) do
                calculate_H2Ag_pes_forces!(n_atoms, coordinates_ang, D)
            end
        end
        
        new(1, potential!, derivative!)
    end
end

function initialize_H2Ag_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call init_pes function
        ccall((:pes_init_, "nn.dylib"),Cvoid ,())
    end
end

function calculate_H2Ag_pes_potential!(n_atoms::Int64, coordinates::Matrix{Float64}, energy::Matrix{Float64})
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:pot0_, "nn.dylib")
    ,Cvoid
    ,(Ref{Int64}
    ,Ref{Float64}
    ,Ptr{Float64}
    )
    , n_atoms, coordinates, energy)
end

function calculate_H2Ag_pes_forces!(n_atoms::Int64, coordinates::Matrix{Float64}, force::Matrix{Float64})
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:dpeshon_, "nn.dylib")
    ,Cvoid
    ,(Ref{Int64}
    ,Ref{Float64}
    ,Ptr{Float64}
    )
    , n_atoms, coordinates, force)
end
