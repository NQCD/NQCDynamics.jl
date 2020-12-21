module EANN_H2_Cu
export EannH2CuModel

using NonadiabaticMolecularDynamics.Atoms
using NonadiabaticMolecularDynamics.Models

"""
Lingjun Zhu, Yaolong Zhang, Liang Zhang, Xueyao Zhou and Bin Jiang
Unified and transferable description of dynamics of H2 dissociative adsorption on multiple copper surfaces via machine learning.
Phys. Chem. Chem. Phys., 2020, 22, 13958--13964
"""

struct EannH2CuModel <: Models.Model

    n_states::UInt
    get_V0::Function
    get_D0::Function

    function EannH2CuModel(path::String, atoms::AtomicParameters)

        function get_V0(R)
            n_atoms = atoms.n_atoms % Int
            res = get_H2Cu_pes_output(path, R, n_atoms, 0, 0, 0)
            return res
        end
        
        function get_D0(R)::Array{Float64} # Vector
            n_atoms = atoms.n_atoms % Int
            res = get_H2Cu_pes_output(path, R, n_atoms, 0, 1, 0)
            return res
        end
        
        new(1, get_V0, get_D0)
    end
end

function initialize_H2Cu_pes()
    # EANN PES: call init_pes function
    ccall((:init_pes_, "nn.dylib"),Cvoid ,())
end

function deallocate_H2Cu_pes() #lib_path::String)
    # EANN PES: call deallocate_all function
    ccall((:deallocate_all_, "nn.dylib"),Cvoid ,())
end

function calculate_h2cu_pes!(coordinates::Array{Float64}, energy::Array{Float64}, force::Array{Float64}, coordinates_type::Int64=0, force_incl::Int64=1)
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
    
    cd(lib_path) do
        dim_atom = 3 # number of coordinates
        n_force=n_atoms-n_constraint # number of forces to calculate (number of unconstrained atoms)
        force=zeros(Float64, dim_atom, n_force) # forces array
        energy=zeros(Float64, 1, 1) # energy value (in 1x1 array because this specific fortran function requires an array not a value)
        
        initialize_H2Cu_pes()

        calculate_h2cu_pes!(coordinates, energy, force, coordinates_type, force_incl)

        deallocate_H2Cu_pes()

        if force_incl == 1
            return force
        else
            return energy
        end
    end
end


end # module
