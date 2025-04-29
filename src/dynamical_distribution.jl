"""
Code uses the multiple dispatch in Julia by defining a certain function multiple times

    So, it depends on your input to automatically choose the appropriate method

    And three dispatches all return the same type DynamicalDistribution

    Therefore, we should have following:

        julia> dist = DynamicalDistribution(10.0, position, (1,1))

        julia> dist.velocity
        NQCDistributions.FixedFill{Float64}(10.0, (1, 1))
"""

"""
    DynamicalDistribution(velocity, position, dims)

Sampleable struct containing distributions for velocity and position.
`dims` determines the size of each sample and should match the size of the system: (ndofs, natoms).

# Example

```jldoctest; setup = :(using Random; Random.seed!(1))
julia> using NQCDistributions: DynamicalDistribution;

julia> d = DynamicalDistribution([[1.0;;], [2.0;;], [3.0;;]], 0.1, (1, 1));

julia> rand(d)
ComponentVector{Float64}(v = [1.0;;], r = [0.1;;])

julia> d[2]
ComponentVector{Float64}(v = [2.0;;], r = [0.1;;])
```
"""
struct DynamicalDistribution{V,R}
    velocity::V
    position::R
    rng::Xoshiro # Random seed generator
    frozen_atoms::Vector{Int} # Indices of atoms which are frozen in place
end

function freeze(velocity_array, frozen_atoms::Vector{Int})
    if isempty(frozen_atoms) # No need to modify the velocity array if no atoms are frozen
        return velocity_array
    else
        for (num, slice) in enumerate(eachslice(velocity_array; dims=2))
            if num in frozen_atoms
                slice .= 0.0 # Set the velocity of the frozen atom to zero
            end
        end
        return velocity_array
    end
end

function DynamicalDistribution(velocity::SampleableComponent, position::SampleableComponent, frozen_atoms::Vector{Int})
    size(velocity) == size(position) || throw(
        DimensionMismatch(
            "`velocity` and `position` sample size does not match: \
            $(size(velocity)) != $(size(position))"
        )
    )
    return DynamicalDistribution(velocity, position, Xoshiro(), frozen_atoms)
end

function DynamicalDistribution(velocity, position, dims::Dims{2}; frozen_atoms=Int[])
    v = SampleableComponent(velocity, dims)
    r = SampleableComponent(position, dims)
    return DynamicalDistribution(v, r, frozen_atoms)
end

function DynamicalDistribution(velocity, position, dims::Dims{3}; classical=Int[], frozen_atoms=Int[])
    v = SampleableComponent(velocity, dims, classical)
    r = SampleableComponent(position, dims, classical)
    return DynamicalDistribution(v, r, frozen_atoms)
end

function Random.rand(rng::AbstractRNG, d::SamplerTrivial{<:DynamicalDistribution})

    if isindexable(d[].velocity) || isindexable(d[].position)
        i = rand(rng, 1:lastindex(d[]))
        return d[][i]
    else
        i = rand(rng, UInt)
        Random.seed!(d[].rng, i)
        return ComponentVector(v=freeze(rand(d[].rng, d[].velocity), d.self.frozen_atoms), r=rand(d[].rng, d[].position))
    end

end

function Base.getindex(d::DynamicalDistribution, i)

    if isindexable(d.velocity) && isindexable(d.position)
        return ComponentVector(v=freeze(d.velocity[i], d.frozen_atoms), r=d.position[i])

    elseif isindexable(d.velocity)
        Random.seed!(d.rng, i)
        return ComponentVector(v=freeze(d.velocity[i], d.frozen_atoms), r=rand(d.rng, d.position))

    elseif isindexable(d.position)
        Random.seed!(d.rng, i)
        return ComponentVector(v=freeze(rand(d.rng, d.velocity), d.frozen_atoms), r=d.position[i])

    else
        Random.seed!(d.rng, i)
        return ComponentVector(v=freeze(rand(d.rng, d.velocity), d.frozen_atoms), r=rand(d.rng, d.position))
    end

end

function Base.lastindex(d::DynamicalDistribution)
    if isindexable(d.velocity) && isindexable(d.position)
        lastindex(d.velocity) == lastindex(d.position) || throw(
            DimensionMismatch("Length of position and velocity distributions do not match.")
        )
        return freeze(lastindex(d.velocity), d.frozen_atoms)

    elseif isindexable(d.velocity)
        return freeze(lastindex(d.velocity), d.frozen_atoms)

    elseif isindexable(d.position)
        return lastindex(d.position)

    else
        return 1

    end
end
Base.length(d::DynamicalDistribution) = lastindex(d)

function Base.iterate(d::DynamicalDistribution, state=1)
    return state > lastindex(d) ? nothing : (d[state], state + 1)
end
