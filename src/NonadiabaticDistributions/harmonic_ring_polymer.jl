using NQCDynamics: RingPolymers
using Distributions: Distributions

struct HarmonicRingPolymer{T,D}
    normal_mode_distribution::D
    normal_mode_transformation::Matrix{T}
end

"""
    HarmonicRingPolymer(ω, β, m, n_beads)

Ring polymer position distribution in a 1D harmonic potential
"""
function HarmonicRingPolymer(ω, β, m, n_beads)
    βₙ = β / n_beads
    ωₙ = 1 / βₙ
    ωₖ = RingPolymers.get_matsubara_frequencies(n_beads, ωₙ)
    σ⁻¹ = sqrt.(βₙ*m .* (ω^2 .+ ωₖ.^2))
    σ = 1 ./ σ⁻¹

    normal_mode_distribution = Distributions.MvNormal(σ)
    normal_mode_transformation = RingPolymers.get_normal_mode_transformation(n_beads)
    HarmonicRingPolymer(normal_mode_distribution, normal_mode_transformation)
end

function Base.rand(rng::AbstractRNG, s::HarmonicRingPolymer)
    R = rand(rng, s.normal_mode_distribution)
    return s.normal_mode_transformation * R
end

function select_item(x::Vector{<:HarmonicRingPolymer{T}}, ::Integer, size::NTuple) where {T}
    length(size) == 3 || throw(error("`size` $size must have 3 dimensions."))
    size[2] == length(x) || throw(error("Sample size does not match distribution."))
    size[3] == length(x[1].normal_mode_distribution) || throw(error("Sample size does not match distribution."))

    out = Array{T,3}(undef, size)
    @views for i in axes(out, 2)
        for j in axes(out, 1)
            out[j,i,:] .= rand(x[i])
        end
    end

    return out
end
