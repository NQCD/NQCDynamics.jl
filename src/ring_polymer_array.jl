using StatsBase: mean
using LinearAlgebra: mul!

export RingPolymerArray
export transform!
export get_centroid

mutable struct RingPolymerArray{T} <: AbstractArray{T,3}
    data::Array{T,3}
    quantum_atoms::Vector{Int}
    classical_atoms::Vector{Int}
    normal::Bool
    function RingPolymerArray{T}(data::Array{T,3}, quantum, classical) where {T}
        new(data, quantum, classical, false)
    end
end

function RingPolymerArray(data::Array{T,3}; quantum=nothing) where {T}
    quantum === nothing && (quantum = collect(axes(data,2)))
    classical = collect(axes(data,2))
    setdiff!(classical, quantum)
    constrain_classical_atoms!(data, classical, 1)
    RingPolymerArray{T}(data, quantum, classical)
end

Base.size(A::RingPolymerArray) = size(A.data)
Base.getindex(A::RingPolymerArray, i::Int, j::Int, k::Int) = A.data[i,j,k]

function Base.setindex!(A::RingPolymerArray, v, i::Int, j::Int, k::Int)
    setindex!(A.data, v, i, j, k)
    j in A.classical_atoms && constrain_classical_atoms!(A, k)
end

function Base.similar(A::RingPolymerArray, ::Type{S}, dims::Dims) where {S}
    if length(dims) == 3
        return RingPolymerArray(Array{S}(undef, dims); quantum=A.quantum_atoms)
    else
        return similar(A.data, S, dims)
    end
end

struct RingPolymerStyle <: Broadcast.AbstractArrayStyle{3} end
Base.BroadcastStyle(::Type{<:RingPolymerArray}) = RingPolymerStyle()
RingPolymerStyle(::Val{3}) = RingPolymerStyle()
RingPolymerStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

Base.similar(bc::Broadcast.Broadcasted{RingPolymerStyle}, ::Type{Eltype}) where {Eltype} =
    RingPolymerArray(similar(Array{Eltype}, axes(bc)))

"""
    constrain_classical_atoms!(A::RingPolymerArray)

Make sure that the classical atoms all have the same values.
"""
function constrain_classical_atoms!(A::Array{T,3}, classical_atoms::AbstractVector, k::Integer) where {T}
    @views for j in classical_atoms
        A[:,j,:] .= A[:,j,k]
    end
end
constrain_classical_atoms!(A::RingPolymerArray, k::Integer) = constrain_classical_atoms!(A.data, A.classical_atoms, k)

"""
Evaluate centroid of ring polymer.
"""
function get_centroid(A::RingPolymerArray{T})::Matrix{T} where {T}
    if A.normal
        centroid = A[:,:,1]
        centroid[:,A.quantum_atoms] ./= sqrt(size(A,3))
    else
        centroid = dropdims(mean(A.data; dims=3); dims=3)
    end
    centroid
end
