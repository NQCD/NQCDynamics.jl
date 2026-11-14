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
    frozen_atoms::AbstractVector # One-hot vector of atoms to freeze velocities of
end

function freeze(velocity_array, frozen_atoms)
    velocity_array_new = velocity_array .* frozen_atoms'
    return velocity_array_new
end

function DynamicalDistribution(velocity::SampleableComponent, position::SampleableComponent, frozen_atoms::AbstractVector)
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
    frozen_onehot = [i ∉ frozen_atoms for i in 1:dims[2]]
    @assert length(frozen_onehot) == dims[2] "Length of frozen_atoms must match number of atoms."
    return DynamicalDistribution(v, r, frozen_onehot)
end

function DynamicalDistribution(velocity, position, dims::Dims{3}; classical=Int[], frozen_atoms=Int[])
    v = SampleableComponent(velocity, dims, classical)
    r = SampleableComponent(position, dims, classical)
    frozen_onehot = [i ∉ frozen_atoms for i in 1:dims[2]]
    @assert length(frozen_onehot) == dims[2] "Length of frozen_atoms must match number of atoms."
    return DynamicalDistribution(v, r, frozen_onehot)
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
        return lastindex(d.velocity)

    elseif isindexable(d.velocity)
        return lastindex(d.velocity)

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
