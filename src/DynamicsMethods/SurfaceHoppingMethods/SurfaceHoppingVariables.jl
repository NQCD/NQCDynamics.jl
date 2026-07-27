using RecursiveArrayTools
using StructArrays
using ArrayInterface

"""
    Due to the need to multiple dispatch this type mirrors that of SurfaceHoppingVariables from the great 
    pacakge RecursiveArrayTools but has been renamed to SurfaceHoppingVariables. This was taken from
    version 3.33.0 so any future changes to SurfaceHoppingVariables will not be mirrored here.
"""
struct SurfaceHoppingVariables{T, A <: RecursiveArrayTools.ArrayPartition{T}, NT <: NamedTuple} <: AbstractVector{T}
    array_partition::A
    names_to_indices::NT
end
SurfaceHoppingVariables(; kwargs...) = SurfaceHoppingVariables(NamedTuple(kwargs))
function SurfaceHoppingVariables(x::NamedTuple)
    names_to_indices = NamedTuple(
        Pair(symbol, index)
            for (index, symbol) in enumerate(keys(x))
    )

    # enforce homogeneity of eltypes
    @assert all(eltype.(values(x)) .== eltype(first(x)))
    T = eltype(first(x))
    S = typeof(values(x))
    return SurfaceHoppingVariables(RecursiveArrayTools.ArrayPartition{T, S}(values(x)), names_to_indices)
end

# Note: overloading `getproperty` means we cannot access `SurfaceHoppingVariables`
# fields except through `getfield` and accessor functions.
RecursiveArrayTools.ArrayPartition(x::SurfaceHoppingVariables) = getfield(x, :array_partition)

function Base.similar(A::SurfaceHoppingVariables)
    return SurfaceHoppingVariables(
        similar(getfield(A, :array_partition)), getfield(A, :names_to_indices)
    )
end

# return SurfaceHoppingVariables when the requested dims still match the partition layout;
# otherwise fall back to the plain backing array of the correct size. RecursiveArrayTools.ArrayPartition's
# own `similar(A, dims)` already does this degradation (it returns a Vector when
# `dims != size(A)`), and we simply propagate that result instead of trying to
# wrap a non-RecursiveArrayTools.ArrayPartition in a SurfaceHoppingVariables (which would hit the inner
# constructor signature `SurfaceHoppingVariables(::A<:RecursiveArrayTools.ArrayPartition, ::NamedTuple)`).
function Base.similar(A::SurfaceHoppingVariables, dims::NTuple{N, Int}) where {N}
    inner = similar(getfield(A, :array_partition), dims)
    inner isa RecursiveArrayTools.ArrayPartition || return inner
    return SurfaceHoppingVariables(inner, getfield(A, :names_to_indices))
end

# similar array partition of common type
@inline function Base.similar(A::SurfaceHoppingVariables, ::Type{T}) where {T}
    return SurfaceHoppingVariables(
        similar(getfield(A, :array_partition), T), getfield(A, :names_to_indices)
    )
end

function Base.similar(A::SurfaceHoppingVariables, ::Type{T}, dims::NTuple{N, Int}) where {T, N}
    inner = similar(getfield(A, :array_partition), T, dims)
    inner isa RecursiveArrayTools.ArrayPartition || return inner
    return SurfaceHoppingVariables(inner, getfield(A, :names_to_indices))
end

# similar array partition with different types
function Base.similar(
        A::SurfaceHoppingVariables, ::Type{T}, ::Type{S}, R::DataType...
    ) where {T, S}
    return SurfaceHoppingVariables(
        similar(getfield(A, :array_partition), T, S, R...), getfield(A, :names_to_indices)
    )
end

Base.Array(x::SurfaceHoppingVariables) = Array(RecursiveArrayTools.ArrayPartition(x))

function Base.zero(x::SurfaceHoppingVariables{T, S, TN}) where {T, S, TN}
    return SurfaceHoppingVariables{T, S, TN}(zero(RecursiveArrayTools.ArrayPartition(x)), getfield(x, :names_to_indices))
end
Base.zero(A::SurfaceHoppingVariables, dims::NTuple{N, Int}) where {N} = zero(A) # ignore dims since named array partitions are vectors

Base.propertynames(x::SurfaceHoppingVariables) = propertynames(getfield(x, :names_to_indices))
function Base.getproperty(x::SurfaceHoppingVariables, s::Symbol)
    return getindex(RecursiveArrayTools.ArrayPartition(x).x, getproperty(getfield(x, :names_to_indices), s))
end

# this enables x.s = some_array.
@inline function Base.setproperty!(x::SurfaceHoppingVariables, s::Symbol, v)
    index = getproperty(getfield(x, :names_to_indices), s)
    return RecursiveArrayTools.ArrayPartition(x).x[index] .= v
end

# print out SurfaceHoppingVariables as a NamedTuple
Base.summary(x::SurfaceHoppingVariables) = string(typeof(x), " with arrays:")
function Base.show(io::IO, m::MIME"text/plain", x::SurfaceHoppingVariables)
    return show(
        io, m, NamedTuple(Pair.(keys(getfield(x, :names_to_indices)), RecursiveArrayTools.ArrayPartition(x).x))
    )
end

Base.size(x::SurfaceHoppingVariables) = size(RecursiveArrayTools.ArrayPartition(x))
Base.length(x::SurfaceHoppingVariables) = length(RecursiveArrayTools.ArrayPartition(x))
# Delegate indexing to the underlying RecursiveArrayTools.ArrayPartition.
# Use concrete index types to avoid invalidating AbstractArray's generic setindex!.
Base.@propagate_inbounds Base.getindex(x::SurfaceHoppingVariables, i::Int) = RecursiveArrayTools.ArrayPartition(x)[i]
Base.@propagate_inbounds Base.setindex!(x::SurfaceHoppingVariables, v, i::Int) = (RecursiveArrayTools.ArrayPartition(x)[i] = v)

# Indexing with non-scalar indices (UnitRange, Vector{Int}, etc.) goes through
# AbstractArray's generic path, which routes via `similar(A, T, dims)`. NAP's
# `similar(::NAP, T, dims)` cannot in general produce a SurfaceHoppingVariables for
# arbitrary `dims` (the partition layout is fixed by `names_to_indices`), so it
# falls back to a plain Vector — making the inferred return type a small Union.
#
# Mirror RecursiveArrayTools.ArrayPartition's `_unsafe_getindex` shortcut at `array_partition.jl:317`:
# allocate the destination directly off the first underlying array and fill it
# via `_unsafe_getindex!`. The result is always a Vector for non-scalar indexing,
# so `x[I]` is type-stable. This matches the v3 indexing semantics (`x[1:end]`
# returns a `Vector`, not a `SurfaceHoppingVariables`); use `similar(x)` /
# `copy(x)` if you want a SurfaceHoppingVariables back.
Base.@propagate_inbounds function Base._unsafe_getindex(
        ::IndexStyle, A::SurfaceHoppingVariables,
        I::Vararg{Union{Real, AbstractArray}, N}
    ) where {N}
    shape = Base.index_shape(I...)
    dest = similar(getfield(A, :array_partition).x[1], shape)
    Base._unsafe_getindex!(dest, A, I...)
    return dest
end
function Base.map(f, x::SurfaceHoppingVariables)
    return SurfaceHoppingVariables(map(f, RecursiveArrayTools.ArrayPartition(x)), getfield(x, :names_to_indices))
end
Base.mapreduce(f, op, x::SurfaceHoppingVariables) = mapreduce(f, op, RecursiveArrayTools.ArrayPartition(x))
# Base.filter(f, x::SurfaceHoppingVariables) = filter(f, RecursiveArrayTools.ArrayPartition(x))

function Base.similar(x::SurfaceHoppingVariables{T, S, NT}) where {T, S, NT}
    return SurfaceHoppingVariables{T, S, NT}(
        similar(RecursiveArrayTools.ArrayPartition(x)), getfield(x, :names_to_indices)
    )
end

# broadcasting
function Base.BroadcastStyle(::Type{<:SurfaceHoppingVariables})
    return Broadcast.ArrayStyle{SurfaceHoppingVariables}()
end
function Base.similar(
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}},
        ::Type{ElType}
    ) where {ElType}
    x = find_SurfaceHoppingVariables(bc)
    return SurfaceHoppingVariables(similar(RecursiveArrayTools.ArrayPartition(x)), getfield(x, :names_to_indices))
end

# when broadcasting with RecursiveArrayTools.ArrayPartition + another array type, the output is the other array type
function Base.BroadcastStyle(
        ::Broadcast.ArrayStyle{SurfaceHoppingVariables}, ::Broadcast.DefaultArrayStyle{1}
    )
    return Broadcast.DefaultArrayStyle{1}()
end

# hook into RecursiveArrayTools.ArrayPartition broadcasting routines
@inline RecursiveArrayTools.npartitions(x::SurfaceHoppingVariables) = RecursiveArrayTools.npartitions(RecursiveArrayTools.ArrayPartition(x))
@inline RecursiveArrayTools.RecursiveArrayTools.unpack(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}},
    i
) = Broadcast.Broadcasted(
    bc.f, RecursiveArrayTools.RecursiveArrayTools.unpack_args(i, bc.args)
)
@inline RecursiveArrayTools.RecursiveArrayTools.unpack(x::SurfaceHoppingVariables, i) = RecursiveArrayTools.unpack(RecursiveArrayTools.ArrayPartition(x), i)

function Base.copy(A::SurfaceHoppingVariables{T, S, NT}) where {T, S, NT}
    return SurfaceHoppingVariables{T, S, NT}(copy(RecursiveArrayTools.ArrayPartition(A)), getfield(A, :names_to_indices))
end

@inline SurfaceHoppingVariables(
    f::F,
    N,
    names_to_indices
) where {
    F <:
    Function,
} = SurfaceHoppingVariables(
    RecursiveArrayTools.ArrayPartition(ntuple(f, Val(N))), names_to_indices
)

@inline function Base.copy(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}})
    N = RecursiveArrayTools.npartitions(bc)
    @inline function f(i)
        return copy(RecursiveArrayTools.unpack(bc, i))
    end
    x = find_SurfaceHoppingVariables(bc)
    return SurfaceHoppingVariables(f, N, getfield(x, :names_to_indices))
end

@inline function Base.copyto!(
        dest::SurfaceHoppingVariables,
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{SurfaceHoppingVariables}}
    )
    N = RecursiveArrayTools.npartitions(dest, bc)
    @inbounds for i in 1:N
        copyto!(getfield(dest, :array_partition).x[i], RecursiveArrayTools.unpack(bc, i))
    end
    return dest
end

#Overwrite ArrayInterface zeromatrix to work with SurfaceHoppingVariabless & implicit solvers within OrdinaryDiffEq
function ArrayInterface.zeromatrix(A::SurfaceHoppingVariables)
    B = RecursiveArrayTools.ArrayPartition(A)
    # Use foldl with explicit init to preserve array type (important for GPU arrays)
    vecs = vec.(B.x)
    rest = Base.tail(vecs)
    x = isempty(rest) ? vecs[1] : foldl(vcat, rest; init = vecs[1])
    return x .* x' .* false
end

# `x = find_SurfaceHoppingVariables(x)` returns the first `SurfaceHoppingVariables` among broadcast arguments.
find_SurfaceHoppingVariables(bc::Base.Broadcast.Broadcasted) = find_SurfaceHoppingVariables(bc.args)
function find_SurfaceHoppingVariables(args::Tuple)
    return find_SurfaceHoppingVariables(find_SurfaceHoppingVariables(args[1]), Base.tail(args))
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

for fn in (:get_velocities, :get_positions, :get_quantum_subsystem)
    @eval function DynamicsUtils.$fn(bc::Base.Broadcast.Broadcasted)
        Base.Broadcast.Broadcasted(bc.f, map(DynamicsUtils.$fn, bc.args))
    end
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
