using StaticArrays
using RecursiveArrayTools

export Momenta
export Positions

struct Momenta{T, S} <: AbstractVectorOfArray{T, 2, Vector{SVector{S, T}}}
    u::Vector{SVector{S, T}}
end
Momenta(vec::Vector{SVector{S, T}}) where {S, T} = Momenta{eltype(T), S}(vec)
Momenta(vec::AbstractVector{<:AbstractVector{T}}) where {T} = Momenta(SVector{length(vec[1]), T}.(vec))
Momenta(mat::AbstractArray{T, 2}) where {T} = Momenta([vec for vec in eachcol(mat)])

struct Positions{T, S} <: AbstractVectorOfArray{T, 2, Vector{SVector{S, T}}}
    u::Vector{SVector{S, T}}
end
Positions(vec::Vector{SVector{S, T}}) where {S, T} = Positions{eltype(T), S}(vec)
Positions(vec::AbstractVector{<:AbstractVector{T}}) where {T} = Positions(SVector{length(vec[1]), T}.(vec))
Positions(mat::AbstractArray{T, 2}) where {T} = Positions([vec for vec in eachcol(mat)])
