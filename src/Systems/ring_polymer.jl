using LinearAlgebra: Symmetric, SymTridiagonal
using RecursiveArrayTools

export RingPolymerArray
export get_bead_masses
export bead_iterator

const RingPolymerArray{T, N} = ArrayPartition{T, NTuple{N, Matrix{T}}} where {T,N}
function RingPolymerArray(A::Vector{Matrix{T}}) where {T}
    RingPolymerArray{T, length(A)}(Tuple(A))
end

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
function RingPolymerParameters{T}(n_beads::Integer, temperature::Real, atom_types::Vector{Symbol}, quantum_nuclei::Vector{Symbol}) where {T}
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

get_bead_masses(n_beads::Integer, masses::Vector) = repeat(masses, inner=n_beads)

bead_iterator(n_beads::Integer, vector::Vector) = Iterators.partition(vector, n_beads)

"""
    get_normal_mode_transformation(n::Int)::Matrix
    
Get the transformation matrix that converts to normal mode coordinates.
"""
function get_normal_mode_transformation(n::Int)::Matrix
    a = ([exp(2π*im/n * j * k) for j=0:n-1, k=0:n-1])
    (real(a) + imag(a)) / sqrt(n)
end