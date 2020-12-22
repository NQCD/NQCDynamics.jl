using LinearAlgebra: Symmetric, SymTridiagonal
using RecursiveArrayTools

export RingPolymerParameters
export get_bead_masses
export bead_iterator
export transform_to_normal_modes!
export transform_from_normal_modes!

struct RingPolymerParameters{T<:AbstractFloat}
    n_beads::UInt
    ω_n::T
    springs::Symmetric{T}
    normal_mode_springs::Vector{T}
    U::Matrix{T}
    quantum_atoms::Vector{UInt}
    function RingPolymerParameters{T}(
        n_beads::Integer, temperature::Real,
        quantum_atoms::Vector{<:Integer}) where {T<:AbstractFloat}

        ω_n = n_beads * temperature
        new(n_beads, ω_n,
            get_spring_matrix(n_beads, ω_n),
            get_normal_mode_springs(n_beads, ω_n), 
            get_normal_mode_transformation(n_beads),
            quantum_atoms)
    end
end

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

"""
    get_normal_mode_transformation(n::Int)::Matrix
    
Get the transformation matrix that converts to normal mode coordinates.
"""
function get_normal_mode_transformation(n::Int)::Matrix
    a = ([exp(2π*im/n * j * k) for j=0:n-1, k=0:n-1])
    (real(a) + imag(a)) / sqrt(n)
end

function transform_to_normal_modes!(p::RingPolymerParameters, R::Array{T,3}, DoF::Integer) where {T}
    @views for i in p.quantum_atoms
        for j=1:DoF
            R[j,i,:] .= p.U'R[j,i,:]
        end
    end
end

function transform_from_normal_modes!(p::RingPolymerParameters, R::Array{T,3}, DoF::Integer) where {T}
    @views for i in p.quantum_atoms
        for j=1:DoF
            R[j,i,:] .= p.U*R[j,i,:]
        end
    end
end

Base.length(beads::RingPolymerParameters) = beads.n_beads
Base.range(beads::RingPolymerParameters) = range(1; length=length(beads))
