using RecursiveArrayTools
using StructArrays


"""
    Due to the need to multiple dispatch this type mirrors that of NamedArrayPartition from the great 
    pacakge RecursiveArrayTools but has been renamed to SurfaceHoppingVariables. This was taken from
    version 3.33.0 so any future changes to NamedArrayPartition will not be mirrored here.
"""
struct SurfaceHoppingVariables{T, A <: ArrayPartition{T}, NT <: NamedTuple} <: AbstractVector{T}
    array_partition::A
    names_to_indices::NT
end

SurfaceHoppingVariables(; kwargs...) = SurfaceHoppingVariables(NamedTuple(kwargs))
function SurfaceHoppingVariables(x::NamedTuple)
    names_to_indices = NamedTuple(Pair(symbol, index)
    for (index, symbol) in enumerate(keys(x)))

    # enforce homogeneity of eltypes
    @assert all(eltype.(values(x)) .== eltype(first(x)))
    T = eltype(first(x))
    S = typeof(values(x))
    return SurfaceHoppingVariables(ArrayPartition{T, S}(values(x)), names_to_indices)
end

# Note: overloading `getproperty` means we cannot access `SurfaceHoppingVariables` 
# fields except through `getfield` and accessor functions.
ArrayPartition(x::SurfaceHoppingVariables) = getfield(x, :array_partition)

function Base.similar(A::SurfaceHoppingVariables)
    SurfaceHoppingVariables(
        similar(getfield(A, :array_partition)), getfield(A, :names_to_indices))
end

# return ArrayPartition when possible, otherwise next best thing of the correct size
function Base.similar(A::SurfaceHoppingVariables, dims::NTuple{N, Int}) where {N}
    SurfaceHoppingVariables(
        similar(getfield(A, :array_partition), dims), getfield(A, :names_to_indices))
end

# similar array partition of common type
@inline function Base.similar(A::SurfaceHoppingVariables, ::Type{T}) where {T}
    SurfaceHoppingVariables(
        similar(getfield(A, :array_partition), T), getfield(A, :names_to_indices))
end

# return ArrayPartition when possible, otherwise next best thing of the correct size
function Base.similar(A::SurfaceHoppingVariables, ::Type{T}, dims::NTuple{N, Int}) where {T, N}
    SurfaceHoppingVariables(
        similar(getfield(A, :array_partition), T, dims), getfield(A, :names_to_indices))
end

# similar array partition with different types
function Base.similar(
        A::SurfaceHoppingVariables, ::Type{T}, ::Type{S}, R::DataType...) where {T, S}
    SurfaceHoppingVariables(
        similar(getfield(A, :array_partition), T, S, R), getfield(A, :names_to_indices))
end

Base.Array(x::SurfaceHoppingVariables) = Array(ArrayPartition(x))

function Base.zero(x::SurfaceHoppingVariables{T, S, TN}) where {T, S, TN}
    SurfaceHoppingVariables{T, S, TN}(zero(ArrayPartition(x)), getfield(x, :names_to_indices))
end
Base.zero(A::SurfaceHoppingVariables, dims::NTuple{N, Int}) where {N} = zero(A) # ignore dims since named array partitions are vectors

Base.propertynames(x::SurfaceHoppingVariables) = propertynames(getfield(x, :names_to_indices))
function Base.getproperty(x::SurfaceHoppingVariables, s::Symbol)
    getindex(ArrayPartition(x).x, getproperty(getfield(x, :names_to_indices), s))
end

# this enables x.s = some_array. 
@inline function Base.setproperty!(x::SurfaceHoppingVariables, s::Symbol, v)
    index = getproperty(getfield(x, :names_to_indices), s)
    ArrayPartition(x).x[index] .= v
end

# print out SurfaceHoppingVariables as a NamedTuple
Base.summary(x::SurfaceHoppingVariables) = string(typeof(x), " with arrays:")
function Base.show(io::IO, m::MIME"text/plain", x::SurfaceHoppingVariables)
    show(
        io, m, NamedTuple(Pair.(keys(getfield(x, :names_to_indices)), ArrayPartition(x).x)))
end

Base.size(x::SurfaceHoppingVariables) = size(ArrayPartition(x))
Base.length(x::SurfaceHoppingVariables) = length(ArrayPartition(x))
Base.getindex(x::SurfaceHoppingVariables, args...) = getindex(ArrayPartition(x), args...)

Base.setindex!(x::SurfaceHoppingVariables, args...) = setindex!(ArrayPartition(x), args...)
function Base.map(f, x::SurfaceHoppingVariables)
    SurfaceHoppingVariables(map(f, ArrayPartition(x)), getfield(x, :names_to_indices))
end
Base.mapreduce(f, op, x::SurfaceHoppingVariables) = mapreduce(f, op, ArrayPartition(x))
# Base.filter(f, x::SurfaceHoppingVariables) = filter(f, ArrayPartition(x))

function Base.similar(x::SurfaceHoppingVariables{T, S, NT}) where {T, S, NT}
    SurfaceHoppingVariables{T, S, NT}(
        similar(ArrayPartition(x)), getfield(x, :names_to_indices))
end

# broadcasting
function Base.BroadcastStyle(::Type{<:SurfaceHoppingVariables})
    Broadcast.ArrayStyle{SurfaceHoppingVariables}()
end
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}},
        ::Type{ElType}) where {ElType}
    x = find_SurfaceHoppingVariables(bc)
    return SurfaceHoppingVariables(similar(ArrayPartition(x)), getfield(x, :names_to_indices))
end

# when broadcasting with ArrayPartition + another array type, the output is the other array tupe
function Base.BroadcastStyle(
        ::Broadcast.ArrayStyle{SurfaceHoppingVariables}, ::Broadcast.DefaultArrayStyle{1})
    Broadcast.DefaultArrayStyle{1}()
end

# hook into ArrayPartition broadcasting routines
@inline RecursiveArrayTools.npartitions(x::SurfaceHoppingVariables) = RecursiveArrayTools.npartitions(ArrayPartition(x))
@inline RecursiveArrayTools.unpack(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}}, i) = Broadcast.Broadcasted(
    bc.f, RecursiveArrayTools.unpack_args(i, bc.args))
@inline RecursiveArrayTools.unpack(x::SurfaceHoppingVariables, i) = RecursiveArrayTools.unpack(ArrayPartition(x), i)

function Base.copy(A::SurfaceHoppingVariables{T, S, NT}) where {T, S, NT}
    SurfaceHoppingVariables{T, S, NT}(deepcopy(ArrayPartition(A)), getfield(A, :names_to_indices))
end

@inline SurfaceHoppingVariables(f::F, N, names_to_indices) where {F <: Function} = SurfaceHoppingVariables(
    ArrayPartition(ntuple(f, Val(N))), names_to_indices)

@inline function Base.copy(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}})
    N = RecursiveArrayTools.npartitions(bc)
    @inline function f(i)
        copy(RecursiveArrayTools.unpack(bc, i))
    end
    x = find_SurfaceHoppingVariables(bc)
    SurfaceHoppingVariables(f, N, getfield(x, :names_to_indices))
end

@inline function Base.copyto!(dest::SurfaceHoppingVariables,
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}})
    N = RecursiveArrayTools.npartitions(dest, bc)
    @inline function f(i)
        copyto!(ArrayPartition(dest).x[i], RecursiveArrayTools.unpack(bc, i))
    end
    ntuple(f, Val(N))
    return dest
end

# `x = find_SurfaceHoppingVariables(x)` returns the first `SurfaceHoppingVariables` among broadcast arguments.
find_SurfaceHoppingVariables(bc::Base.Broadcast.Broadcasted) = find_SurfaceHoppingVariables(bc.args)
function find_SurfaceHoppingVariables(args::Tuple)
    find_SurfaceHoppingVariables(find_SurfaceHoppingVariables(args[1]), Base.tail(args))
end
find_SurfaceHoppingVariables(x) = x
find_SurfaceHoppingVariables(::Tuple{}) = nothing
find_SurfaceHoppingVariables(x::SurfaceHoppingVariables, rest) = x
find_SurfaceHoppingVariables(::Any, rest) = find_SurfaceHoppingVariables(rest)


DynamicsUtils.get_velocities(u::SurfaceHoppingVariables) = u.v
DynamicsUtils.get_positions(u::SurfaceHoppingVariables) = u.r
function DynamicsUtils.get_quantum_subsystem(u::SurfaceHoppingVariables)
    T = eltype(u)
    real::Matrix{T} = u.σreal
    imag::Matrix{T} = u.σimag
    return StructArray{Complex{T}}((real, imag))
end

function LinearAlgebra.mul!(C::SurfaceHoppingVariables, A, v::SurfaceHoppingVariables)
    σ = DynamicsUtils.get_quantum_subsystem(v)
    C = DynamicsUtils.get_quantum_subsystem(C)
    return mul!(C, A, σ)
end

OrdinaryDiffEq.OrdinaryDiffEqCore._vec(x::SurfaceHoppingVariables) = x.σreal[:,1]

import Base.(*)

function (*)(A::AbstractMatrix{T}, x::SurfaceHoppingVariables) where {T}
    σ = DynamicsUtils.get_quantum_subsystem(x)
    mul!(similar(σ), A, σ)
end

function (*)(x::SurfaceHoppingVariables, A::AbstractMatrix{T}, ) where {T}
    σ = DynamicsUtils.get_quantum_subsystem(x)
    mul!(similar(σ), A, σ)
end

@inline function Base.materialize!(dest::SurfaceHoppingVariables, bc::Base.Broadcast.Broadcasted{<:Any})
    SurfaceHoppingVariables_broadcast(dest, bc, bc.style)
end

@inline function SurfaceHoppingVariables_broadcast(dest, bc, ::Any)
    return Base.materialize!(Base.Broadcast.combine_styles(dest, bc), dest, bc)
end

@inline function SurfaceHoppingVariables_broadcast(dest, bc, ::Base.Broadcast.ArrayStyle{SurfaceHoppingVariables})
    new_des = DynamicsUtils.get_quantum_subsystem(dest)
    original = DynamicsUtils.get_quantum_subsystem(bc.args[1])
    new_des .= original
end

@inline function SurfaceHoppingVariables_broadcast(dest, bc, ::StructArrays.StructArrayStyle)
    new_des = DynamicsUtils.get_quantum_subsystem(dest)
    original = bc.args[1]
    new_des .= original
end

# -------------------------------------------- TESTING ------------------------------------------- #
#= function int(du,u,p,t)
    du.y .= u.y .+ t/10
    du.x .= u.x.+1.0
end

function int2(du,u,p,t)
    du.x[2] .= u.x[2] .+ t/10
    du.x[1] .= u.x[1].+1.0
end

prob = ODEProblem(int, u0, (0.0,2.0))
sol = solve(prob, Tsit5()) =#
# ------------------------------------------------------------------------------------------------ #
