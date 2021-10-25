using Distributions: Distributions
using UnitfulAtomic: austrip
using Random: Random

struct BoltzmannVelocityDistribution{T,S} <: VelocityDistribution
    dist::Distributions.MvNormal{T}
    size::NTuple{S,Int}
end

"""
    BoltzmannVelocityDistribution(temperature, masses)

Obtain velocities sampled from the Boltzmann distribution for given temperature and masses.

This is effectively just a multivariate normal distribution but when combined with the
[`DynamicalDistribution`](@ref) it ensure samples are generated with the correct size.
If using for ring polymer dynamics, ensure the temperature is given as nbeads*temperature.
"""
function BoltzmannVelocityDistribution(temperature, masses, size)
    dist = Distributions.MvNormal(sqrt.(austrip(temperature) ./ masses))
    BoltzmannVelocityDistribution(dist, size)
end

function Base.rand(rng::AbstractRNG, s::BoltzmannVelocityDistribution{T,S}) where {T,S}
    out = Array{T,S}(undef, s.size...)
    if S == 2
        @views for i=1:size(out, 1)
            Random.rand!(rng, s.dist, out[i,:])
        end
    elseif S == 3
        @views for i=1:size(out, 3)
            for j=1:size(out, 1)
                Random.rand!(rng, s.dist, out[j,:,i])
            end
        end
    end
    return out
end

select_item(x::BoltzmannVelocityDistribution, ::Integer, ::NTuple) = rand(x)
