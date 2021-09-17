
using Random: AbstractRNG, rand!
using Distributions: Distributions, Continuous, Sampleable, Univariate
using HDF5: h5open

struct DynamicalVariate <: Distributions.VariateForm end

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
Base.eltype(s::DynamicalDistribution{<:AbstractArray,R} where {R}) = eltype(s.velocity[1])
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
select_item(x::Vector{<:AbstractArray}, i::Integer, ::NTuple) = x[i]
select_item(x::Vector{<:Number}, i::Integer, size::NTuple) = fill(x[i], size)
function select_item(x::AbstractString, i::Integer, ::NTuple)
    h5open(x, "r") do fid
        return read(fid, string(i))
    end
end

# Sampled selection
select_item(x::Sampleable{Univariate}, ::Integer, size::NTuple) = rand(x, size)
function select_item(x::Vector{<:Sampleable{Univariate}}, ::Integer, size::NTuple)
    if length(x) == size[2]
        if length(size) == 3
            tmpsize = (size[1], 1, size[3])
        elseif length(size) == 2
            tmpsize = (size[1], 1)
        end
        a = rand.(x, Ref(tmpsize))
        return cat(a...; dims=2)
    else
        return throw(ErrorException(
            "Distribution size does not match sample size. natoms != length(distribution)"
            ))
    end
end

# Deterministic selections
select_item(x::Real, ::Integer, size::NTuple) = fill(x, size)
select_item(x::Matrix, ::Integer, ::NTuple) = x
select_item(x::AbstractArray{T,3}, ::Integer, ::NTuple) where T = x

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
