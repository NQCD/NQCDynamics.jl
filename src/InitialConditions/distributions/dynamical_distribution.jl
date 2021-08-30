
export DynamicalDistribution

using Random: AbstractRNG, rand!
using Distributions
using UnitfulAtomic
using HDF5

struct DynamicalVariate <: VariateForm end

struct DynamicalDistribution{V,R,S,N} <: Sampleable{DynamicalVariate,Continuous}
    velocity::V
    position::R
    size::NTuple{S,Int}
    state::N
    type::Symbol
end
DynamicalDistribution(velocity, position, size; state=0, type=:adiabatic) =
    DynamicalDistribution(velocity, position, size, state, type)

Base.eltype(s::DynamicalDistribution{<:Sampleable,R}) where {R} = eltype(s.velocity)
Base.eltype(s::DynamicalDistribution{<:AbstractArray,R} where {R}) = eltype(austrip.(s.velocity[1]))
Base.size(s::DynamicalDistribution) = s.size

function Distributions.rand(rng::AbstractRNG, s::Sampleable{DynamicalVariate})
    Distributions._rand!(rng, s, [Array{eltype(s)}(undef, size(s)) for i=1:2])
end

function maxindex(s::DynamicalDistribution)
    if s.velocity isa Vector
        return length(s.velocity)
    elseif s.position isa Vector
        return length(s.position)
    elseif s.velocity isa String
        h5open(s.velocity, "r") do fid
            return length(fid)
        end
    elseif s.position isa String
        h5open(s.position, "r") do fid
            return length(fid)
        end
    else
        return 1
    end
end

function Distributions._rand!(rng::AbstractRNG, s::DynamicalDistribution, x::Vector{<:Array})
    i = rand(rng, 1:maxindex(s))
    x[1] .= select_item(s.velocity, i, s.size)
    x[2] .= select_item(s.position, i, s.size)
    x
end

pick(s::DynamicalDistribution, i::Integer) = [select_item(s.velocity, i, s.size), select_item(s.position, i, s.size)]

# Indexed selections
select_item(x::Vector{<:AbstractArray}, i::Integer, ::NTuple) = austrip.(x[i])
select_item(x::Vector{<:Number}, i::Integer, size::NTuple) = fill(austrip.(x[i]), size)
function select_item(x::AbstractString, i::Integer, ::NTuple)
    h5open(x, "r") do fid
        return read(fid, string(i))
    end
end

# Sampled selection
select_item(x::Sampleable{Univariate}, ::Integer, size::NTuple) = austrip.(rand(x, size))

# Deterministic selections
select_item(x::Real, ::Integer, size::NTuple) = austrip.(fill(x, size))
select_item(x::Matrix, ::Integer, ::NTuple) = austrip.(x)
select_item(x::AbstractArray{T,3}, ::Integer, ::NTuple) where T = austrip.(x)

function Base.show(io::IO, s::DynamicalDistribution) 
    print(io, "DynamicalDistribution with size: ", size(s))
end

function Base.show(io::IO, ::MIME"text/plain", s::DynamicalDistribution)
    print(io, "DynamicalDistribution:\n  ",
          "size: ", size(s), "\n  ",
          "state: ", s.state, "\n  ",
          "type: ", s.type)
end

function write_hdf5(filename, data)
    h5open(filename, "w") do fid
        for (i,d) in enumerate(data)
            fid[string(i)] = d
        end
    end
end
