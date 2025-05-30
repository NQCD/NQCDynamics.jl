using LinearAlgebra: LinearAlgebra, Symmetric, SymTridiagonal, I

using NQCBase: NQCBase, PeriodicCell, InfiniteCell
using RingPolymerArrays: RingPolymerArrays, NormalModeTransformation

struct RingPolymerParameters{T<:AbstractFloat}
    n_beads::Int
    ω_n::T
    ω_n²::T
    springs::Symmetric{T,Matrix{T}}
    normal_mode_springs::Vector{T}
    transformation::NormalModeTransformation{T}
    quantum_atoms::Vector{Int}
    function RingPolymerParameters{T}(
        n_beads::Integer, temperature::Real,
        quantum_atoms::Vector{<:Integer}) where {T<:AbstractFloat}

        ω_n = n_beads * temperature
        new(n_beads, ω_n, ω_n^2,
            get_spring_matrix(n_beads, ω_n),
            get_normal_mode_springs(n_beads, ω_n), 
            NormalModeTransformation{T}(n_beads),
            quantum_atoms,
        )
    end
end

Base.broadcastable(beads::RingPolymerParameters) = Ref(beads)

"""Constructor for choosing specific elements to be quantum."""
function RingPolymerParameters{T}(n_beads::Integer, temperature::Real, atom_types::AbstractVector{Symbol}, quantum_nuclei::Vector{Symbol}) where {T}
    quantum_atoms = findall(in(quantum_nuclei), atom_types)
    RingPolymerParameters{T}(n_beads, temperature, quantum_atoms)
end

"""Constructor for the case where all nuclei are quantum."""
function RingPolymerParameters{T}(n_beads::Integer, temperature::Real, n_atoms::Integer) where {T}
    RingPolymerParameters{T}(n_beads, temperature, collect(1:n_atoms))
end

"""
    get_L(n_beads, mass, ω_n)

Return the Circulant symmetric matrix for the ring polymer springs.
"""
function get_spring_matrix(n_beads::Integer, ω_n::Real)::Symmetric
    if n_beads == 1
        spring = zeros(1, 1)
    elseif n_beads == 2
        spring = [2 -2; -2 2]
    else
        spring = SymTridiagonal(fill(2, n_beads), fill(-1, n_beads-1))
        spring = convert(Matrix, spring)
        spring[end,1] = spring[1, end] = -1
    end
    Symmetric(spring .*  ω_n^2 / 2)
end

get_normal_mode_springs(n_beads::Integer, ω_n::Real) = get_matsubara_frequencies(n_beads, ω_n) .^2 / 2
get_matsubara_frequencies(n::Integer, ω_n::Real) = 2ω_n*sin.((0:n-1)*π/n)

nbeads(beads::RingPolymerParameters) = beads.n_beads
Base.length(beads::RingPolymerParameters) = nbeads(beads)
Base.range(beads::RingPolymerParameters) = range(1; length=length(beads))

"""
    cayley_propagator(beads::RingPolymerParameters{T}, dt::Real; half::Bool=true) where {T}

J. Chem. Phys. 151, 124103 (2019); doi: 10.1063/1.5120282
"""
function cayley_propagator(beads::RingPolymerParameters{T}, dt::Real; half::Bool=true) where {T}

    cay(dtA::Matrix)::Matrix = LinearAlgebra.inv(I - dtA/2) * (I + dtA/2)

    ω_k = get_matsubara_frequencies(length(beads), beads.ω_n)
    prop = [Array{T}(undef, 2, 2) for i=1:length(beads)]
    for (i, ω) in enumerate(ω_k)
        A = [0 1; -ω^2 0]
        prop[i] .= half ? real.(sqrt(cay(dt.*A))) : cay(dt.*A)
    end
    prop
end

"""
    get_spring_energy(beads::RingPolymerParameters, masses, R)
    
Calculate the ring polymer spring potential.
"""
function get_spring_energy(beads::RingPolymerParameters, masses, R)
    E = zero(eltype(R))

    for i in beads.quantum_atoms
        for j in axes(R, 1)
            E += masses[i] * (R[j, i, end] - R[j, i, 1])^2
        end
    end

    for bead=1:nbeads(beads)-1
        for i in beads.quantum_atoms # Only for quantum nuclei
            for j in axes(R, 1)
                E += masses[i] * (R[j, i, bead] - R[j, i, bead+1])^2
            end
        end
    end

    return E * beads.ω_n^2/2
end

function NQCBase.apply_cell_boundaries!(cell::PeriodicCell, R::AbstractArray{T,3}, beads::RingPolymerParameters) where {T}
    RingPolymerArrays.transform_to_normal_modes!(R, beads.transformation)
    R[:,beads.quantum_atoms,1] ./= sqrt(length(beads))
    @views NQCBase.apply_cell_boundaries!(cell, R[:,:,1])
    R[:,beads.quantum_atoms,1] .*= sqrt(length(beads))
    RingPolymerArrays.transform_from_normal_modes!(R, beads.transformation)
end

NQCBase.apply_cell_boundaries!(::InfiniteCell, ::AbstractArray{T,3}, ::RingPolymerParameters) where {T} = nothing
