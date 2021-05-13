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
    function RingPolymerArray{T}(data::Array{T,3}, quantum, classical, normal) where {T}
        constrain_classical_atoms!(data, classical)
        new(data, quantum, classical, normal)
    end
end

function RingPolymerArray(data::Array{T,3}; quantum=nothing, normal=false) where {T}
    quantum === nothing && (quantum = collect(axes(data,2)))
    classical = collect(axes(data,2))
    setdiff!(classical, quantum)
    RingPolymerArray{T}(data, quantum, classical, normal)
end

Base.size(A::RingPolymerArray) = size(A.data)
Base.getindex(A::RingPolymerArray, i::Int, j::Int, k::Int) = A.data[i,j,k]

function Base.setindex!(A::RingPolymerArray, v, i::Int, j::Int, k::Int)
    if j in A.classical_atoms
        if k == 1
            setindex!(A.data, v, i, j, k)
            set_classical!(A, i, j)
        end
    else
        setindex!(A.data, v, i, j, k)
    end
end

function Base.similar(A::RingPolymerArray, ::Type{S}, dims::Dims) where {S}
    if length(dims) == 3
        return RingPolymerArray{S}(Array{S}(undef, dims), A.quantum_atoms, A.classical_atoms, A.normal)
    else
        return similar(A.data, S, dims)
    end
end

struct RingPolymerStyle <: Broadcast.AbstractArrayStyle{3} end
Base.BroadcastStyle(::Type{<:RingPolymerArray}) = RingPolymerStyle()
RingPolymerStyle(::Val{3}) = RingPolymerStyle()
RingPolymerStyle(::Val{N}) where {N} = Broadcast.DefaultArrayStyle{N}()

function Base.similar(bc::Broadcast.Broadcasted{RingPolymerStyle}, ::Type{Eltype}) where {Eltype}
    A = find_rpa(bc)
    if A === nothing
        return similar(Array{Eltype}, axes(bc))
    else
        return similar(A)
    end
end

"`A = find_rpa(As)` return the first RingPolymerArray among the arguments."
find_rpa(bc::Base.Broadcast.Broadcasted) = find_rpa(bc.args)
find_rpa(args::Tuple) = find_rpa(find_rpa(args[1]), Base.tail(args))
find_rpa(x) = x
find_rpa(::Tuple{}) = nothing
find_rpa(a::RingPolymerArray, rest) = a
find_rpa(::Any, rest) = find_rpa(rest)

"""
    constrain_classical_atoms!(A::AbstractArray{T,3}, classical_atoms) where {T}

Make sure that the classical atoms all have the same values.
"""
function constrain_classical_atoms!(A::AbstractArray{T,3}, classical_atoms) where {T}
    for j in classical_atoms
        for i in axes(A,1)
            set_classical!(A, i, j)
        end
    end
end

"Copy data from first index into replicas. All classical atoms have repeated values."
set_classical!(A::AbstractArray{T,3}, i, j) where {T} = A[i,j,2:end] .= A[i,j,1]
set_classical!(A::RingPolymerArray, i, j) = set_classical!(A.data, i, j)

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
