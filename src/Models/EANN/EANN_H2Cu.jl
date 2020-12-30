export EANN_H₂Cu

"""
Lingjun Zhu, Yaolong Zhang, Liang Zhang, Xueyao Zhou and Bin Jiang
Unified and transferable description of dynamics of H2 dissociative adsorption on multiple copper surfaces via machine learning.
Phys. Chem. Chem. Phys., 2020, 22, 13958--13964
"""

struct EANN_H₂Cu <: AdiabaticModel

    n_states::UInt
    potential!::Function
    derivative!::Function

    function EANN_H₂Cu(path::String, atoms::Atoms)

        initialize_H2Cu_pes(path)

        n_atoms = convert(Int, length(atoms))

        function potential!(V::AbstractMatrix, R::AbstractMatrix)
            V .= get_H2Cu_pes_output(path, R, n_atoms, 0, 0, 0)
        end
        
        function derivative!(D::AbstractMatrix, R::AbstractMatrix)
            D .= get_H2Cu_pes_output(path, R, n_atoms, 0, 1, 0)
        end
        
        new(1, potential!, derivative!)
    end
end

function initialize_H2Cu_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call init_pes function
        ccall((:init_pes_, "nn.dylib"),Cvoid ,())
    end
end

function deallocate_H2Cu_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call deallocate_all function
        ccall((:deallocate_all_, "nn.dylib"),Cvoid ,())
    end
end

function calculate_H2Cu_pes!(coordinates::Array{Float64}, energy::Array{Float64}, force::Array{Float64}, coordinates_type::Int64=0, force_incl::Int64=1)
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:eann_out_, "nn.dylib")
    ,Cvoid 
    ,(Ref{Int64}
    ,Ref{Int64}
    ,Ref{Float64}
    ,Ptr{Float64}
    ,Ptr{Float64} 
    )
    , coordinates_type, force_incl, coordinates, energy, force)
end

function get_H2Cu_pes_output(lib_path::String, coordinates::Array{Float64}, n_atoms::Int64, coordinates_type::Int64=0, force_incl::Int64=1, n_constraint::Int64=0)
    
    # have to use cd(path) because of the issue with calling the ccall functions using an entire path and issues with @eval
    cd(lib_path) do
        dim_atom = 3 # number of coordinates
        n_force=n_atoms-n_constraint # number of forces to calculate (number of unconstrained atoms)
        force=zeros(Float64, dim_atom, n_force) # forces array
        energy=zeros(Float64, 1, 1) # energy value (in 1x1 array because this specific fortran function requires an array not a value)
        coordinates_ang = copy(ustrip(auconvert.(u"Å", coordinates))) # converts coordinates into Angstrom

        calculate_H2Cu_pes!(coordinates_ang, energy, force, coordinates_type, force_incl)

        if force_incl == 1
            return force
        else
            return energy
        end
    end
end
