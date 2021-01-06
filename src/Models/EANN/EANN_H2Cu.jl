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

        function potential!(V::AbstractVector, R::AbstractMatrix)
            coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
            cd(path) do
                calculate_H2Cu_pes!(coordinates_ang, reshape(V, (1,1)), zero(R), 0, 0)
            end
        end
        
        function derivative!(D::AbstractMatrix, R::AbstractMatrix)
            coordinates_ang = copy(ustrip(auconvert.(u"Å", R)))
            cd(path) do
                calculate_H2Cu_pes!(coordinates_ang, zeros(1,1), D, 0, 1)
            end
            D .*= -1
        end
        
        new(1, potential!, derivative!)
    end
end

function initialize_H2Cu_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call init_pes function
        ccall((:init_pes_, "h2cu_pes.dylib"), Cvoid, ())
    end
end

function deallocate_H2Cu_pes(lib_path::String)
    cd(lib_path) do
        # EANN PES: call deallocate_all function
        ccall((:deallocate_all_, "h2cu_pes.dylib"), Cvoid, ())
    end
end

function calculate_H2Cu_pes!(coordinates::Array{Float64}, energy::Array{Float64}, force::Array{Float64}, coordinates_type::Int64=0, force_incl::Int64=1)
    # EANN PES: get energy and force from EANN to vars: energy and force
    ccall((:eann_out_, "h2cu_pes.dylib"),
        Cvoid,
        (Ref{Int64}, Ref{Int64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}),
        coordinates_type, force_incl, coordinates, energy, force)
end
