
using Random: AbstractRNG
using Distributions: Univariate, Sampleable
using HDF5: h5open
using ComponentArrays: ComponentVector

"""
    DynamicalDistribution(velocity, position, size)

Sampleable that provides positions and velocities sampled from a variety of distributions.

# Arguments

**`velocity`** and **`position`** provide the velocities and positions and can be any type for
which `select_item` has been implemented.

**`size`** denotes to the size of the samples that should be produced.

# Keywords

# Example

```jldoctest; setup = :(using Random; Random.seed!(1))
using NonadiabaticMolecularDynamics.InitialConditions
using Distributions

d = DynamicalDistribution(5.0, Normal(), (1, 1))
rand(d)

# output

ComponentVector{Float64}(v = [5.0], r = [-0.3170409357632898])
```
"""
struct DynamicalDistribution{V,R,S} <: NuclearDistribution
    velocity::V
    position::R
    size::NTuple{S,Int}
end

Base.size(s::DynamicalDistribution) = s.size
Base.length(s::DynamicalDistribution) = maxindex(s)

function Base.rand(rng::AbstractRNG, s::DynamicalDistribution)
    i = rand(rng, 1:maxindex(s))
    return pick(s, i)
end

function maxindex(s::DynamicalDistribution)::Int
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
        return 1 end
end

pick(s::DynamicalDistribution, i::Integer) = ComponentVector(v=select_item(s.velocity, i, s.size), r=select_item(s.position, i, s.size))

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
    print(io, "DynamicalDistribution:\n  ", "size: ", size(s))
end

function write_hdf5(filename, data)
    h5open(filename, "w") do fid
        for (i,d) in enumerate(data)
            fid[string(i)] = d
        end
    end
end

function Base.iterate(s::DynamicalDistribution, state=1)
    state > maxindex(s) ? nothing : (pick(s, state), state+1)
end
