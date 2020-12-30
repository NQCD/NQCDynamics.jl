
module EANN_H2Ag

export EannH2AgModel

using NonadiabaticMolecularDynamics.Atoms
using NonadiabaticMolecularDynamics.Models

"""
J. Phys. Chem. Lett. 2019, 10, 4962−4967
J. Phys. Chem. C 2020, 124, 186−195
"""

struct EannH2AgModel <: Models.Model

    n_states::UInt
    get_V0::Function
    get_D0::Function

    function EannH2AgModel(path::String, atoms::AtomicParameters)

        function get_V0(R)
            n_atoms = atoms.n_atoms % Int
            res = get_H2Ag_pes_output(path, R, n_atoms, 0)
            return res
        end
        
        function get_D0(R)::Array{Float64} # Vector
            n_atoms = atoms.n_atoms % Int
            res = get_H2Ag_pes_output(path, R, n_atoms, 1)
            return res
        end
        
        new(1, get_V0, get_D0)
    end
end

function initialize_H2Ag_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call init_pes function
        ccall((:pes_init_, "nn.dylib"),Cvoid ,())
    end
end

function calculate_H2Ag_pes_potential!(n_atoms::Int64, coordinates::Array{Float64}, energy::Array{Float64})
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:pot0_, "nn.dylib")
    ,Cvoid
    ,(Ref{Int64}
    ,Ref{Float64}
    ,Ptr{Float64}
    )
    , n_atoms, coordinates, energy)
end

function calculate_H2Ag_pes_forces!(n_atoms::Int64, coordinates::Array{Float64}, force::Array{Float64})
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:dpeshon_, "nn.dylib")
    ,Cvoid
    ,(Ref{Int64}
    ,Ref{Float64}
    ,Ptr{Float64}
    )
    , n_atoms, coordinates, force)
end

function get_H2Ag_pes_output(lib_path::String, coordinates::Array{Float64}, n_atoms::Int64, force_incl::Int64=1)
    
    cd(lib_path) do
        dim_atom = 3 # number of coordinates
        n_force=n_atoms #-n_constraint # number of forces to calculate (number of unconstrained atoms)
        force=zeros(Float64, dim_atom, n_force) # forces array
        energy=zeros(Float64, 1, 1) # energy value (in 1x1 array because this specific fortran function requires an array not a value)
        coordinates_ang = copy(ustrip(auconvert.(u"Å", coordinates))) # converts coordinates into Angstrom

        if force_incl == 1
            calculate_H2Ag_pes_forces!(n_atoms, coordinates_ang, force)
            return force
        else
            calculate_H2Ag_pes_potential!(n_atoms, coordinates_ang, energy)
            return energy
        end
    end
end


end # module
