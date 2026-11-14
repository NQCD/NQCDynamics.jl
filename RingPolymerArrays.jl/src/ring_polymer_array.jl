
"""
    IndexMapping

Small struct containing the index mappings from the original indexing into the
internal `RingPolymerArray` storage format.
"""
struct IndexMapping
    classical::OrderedDict{Int,Int}
    quantum::OrderedDict{Int,Int}
end

function IndexMapping(quantum::Vector{Int}, classical::Vector{Int})
    classical_map = OrderedDict(classical .=> 1:length(classical))
    quantum_map = OrderedDict(quantum .=> 1:length(quantum))
    return IndexMapping(classical_map, quantum_map)
end

"""
    RingPolymerArray{T} <: AbstractArray{T,3}

Array for representing ring polymer systems.

The system is partioned into `classical_atoms` and `quantum_atoms`,
where the `classical_atoms` have only a single bead and the `quantum_atoms` have many.
The total number of beads is equal to the size of the third dimension.
"""
struct RingPolymerArray{T} <: AbstractArray{T,3}
    classical_atoms::Matrix{T}
    quantum_atoms::Vector{Matrix{T}}
    index_map::IndexMapping
    dims::Dims{3}
end

function RingPolymerArray{T}(A::AbstractArray; classical=Int[]) where {T}
    return RingPolymerArray(T.(A); classical)
end

function RingPolymerArray(A::AbstractArray{T,3}; classical=Int[]) where {T}
    B = RingPolymerArray{T}(undef, size(A); classical)
    copy!(B, A)
    return B
end

function RingPolymerArray(A::AbstractMatrix{T}; nbeads, classical=Int[]) where {T}
    B = RingPolymerArray{T}(undef, (size(A)..., nbeads); classical)
    for bead in eachbead(B)
        copy!(bead, A)
    end
    return B
end

function RingPolymerArray{T}(::UndefInitializer, dims::Dims{3}; classical::AbstractVector{<:Integer}) where {T}
    quantum = find_quantum_indices(dims[2], classical)
    index_map = IndexMapping(quantum, classical)
    classical_atoms = Matrix{T}(undef, dims[1], length(classical))
    quantum_atoms = [Matrix{T}(undef, dims[1], length(quantum)) for _ in 1:dims[3]]
    return RingPolymerArray{T}(classical_atoms, quantum_atoms, index_map, dims)
end

"""
    find_quantum_indices(natoms::Integer, classical::AbstractVector{<:Integer})

Returns all indices not included in `classical`.
"""
function find_quantum_indices(natoms::Integer, classical::AbstractVector{<:Integer})
    quantum = collect(1:natoms)
    setdiff!(quantum, classical)
    return quantum
end

function Base.similar(A::RingPolymerArray, ::Type{T}, dims::Dims{3}) where {T}
    return RingPolymerArray{T}(undef, dims; classical=collect(classicalindices(A)))
end

function Base.similar(::RingPolymerArray, ::Type{T}, dims::Dims) where {T}
    return Array{T}(undef, dims)
end

Base.size(A::RingPolymerArray) = A.dims

function Base.getindex(A::RingPolymerArray, i, j, k)
    j in classicalindices(A) ? get_classical_index(A, i, j, k) : get_quantum_index(A, i, j, k)
end

get_classical_index(A::RingPolymerArray, i, j, k) = A.classical_atoms[i, A.index_map.classical[j]]
get_quantum_index(A::RingPolymerArray, i, j, k) = A.quantum_atoms[k][i, A.index_map.quantum[j]]

function Base.setindex!(A::RingPolymerArray, v, i, j, k)
    j in classicalindices(A) ? set_classical_index!(A, v, i, j, k) : set_quantum_index!(A, v, i, j, k)
end

function set_classical_index!(A::RingPolymerArray, v, i, j, k)
    k == 1 && setindex!(A.classical_atoms, v, i, A.index_map.classical[j])
end

function set_quantum_index!(A::RingPolymerArray, v, i, j, k)
    setindex!(A.quantum_atoms[k], v, i, A.index_map.quantum[j])
end

quantumindices(A::RingPolymerArray) = keys(A.index_map.quantum)
classicalindices(A::RingPolymerArray) = keys(A.index_map.classical)

"""
    eachbead(A::AbstractArray{T,3}) where {T}

Iterate views of each bead.
Slices the array along the first two (dofs, atoms) dimensions.
"""
eachbead(A::AbstractArray{T,3}) where {T} = (view(A, :, :, i) for i in axes(A, 3))

"""
    eachdof(A::AbstractArray{T,3}) where {T}

Iterate over every degree of freedom for all beads.
Slices the array along the third (bead) dimension.
"""
eachdof(A::AbstractArray{T,3}) where {T} = (view(A, i, j, :) for i in axes(A, 1), j in axes(A, 2))

"""
    get_centroid(A::AbstractArray{T,3}) where {T}

Get the ring polymer centroid by averaging bead coordinates.

This assumes that the array is not in normal mode coordinates.
"""
function get_centroid(A::AbstractArray{T,3}) where {T}
    centroid = zeros(T, size(A,1), size(A,2))
    get_centroid!(centroid, A)
    return centroid
end

function get_centroid!(centroid::Matrix{T}, A::AbstractArray{T,3}) where {T}
    fill!(centroid, zero(eltype(A)))
    @inbounds for i in axes(A,3)
        for j in axes(A,2)
            for k in axes(A,1)
                centroid[k,j] += A[k,j,i]
            end
        end
    end
    @inbounds for I in eachindex(centroid)
        centroid[I] /= size(A,3)
    end
    return centroid
end

struct RingPolymerStyle <: Broadcast.AbstractArrayStyle{3} end
RingPolymerStyle(::Val{3}) = RingPolymerStyle()
Base.BroadcastStyle(::Type{<:RingPolymerArray}) = RingPolymerStyle()

function Base.similar(bc::Broadcast.Broadcasted{RingPolymerStyle}, ::Type{ElType}) where {ElType}
    A = find_aac(bc)
    return similar(A)
end

# https://docs.julialang.org/en/v1/manual/interfaces/#man-interfaces-broadcasting
find_aac(bc::Base.Broadcast.Broadcasted) = find_aac(bc.args)
find_aac(args::Tuple) = find_aac(find_aac(args[1]), Base.tail(args))
find_aac(x) = x
find_aac(::Tuple{}) = nothing
find_aac(a::RingPolymerArray, rest) = a
find_aac(::Any, rest) = find_aac(rest)
