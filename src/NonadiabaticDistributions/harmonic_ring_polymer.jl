using NQCDynamics: RingPolymers
using Distributions: Distributions
using RingPolymerArrays: RingPolymerArrays, NormalModeTransformation
using LinearAlgebra: mul!

struct HarmonicRingPolymer{T,D}
    normal_mode_distribution::D
    normal_mode_transformation::NormalModeTransformation{T}
end

"""
    HarmonicRingPolymer(ω, β, m, n_beads)

Ring polymer position distribution in a 1D harmonic potential
"""
function HarmonicRingPolymer(ω, β, m, n_beads; centre=0)
    βₙ = β / n_beads
    ωₙ = 1 / βₙ
    ωₖ = RingPolymers.get_matsubara_frequencies(n_beads, ωₙ)
    σ⁻¹ = sqrt.(βₙ*m .* (ω^2 .+ ωₖ.^2))
    σ = 1 ./ σ⁻¹

    μ = [centre * sqrt(n_beads); zeros(n_beads-1)]
    normal_mode_distribution = Distributions.MvNormal(μ, σ)
    normal_mode_transformation = RingPolymerArrays.NormalModeTransformation{Float64}(n_beads)
    HarmonicRingPolymer(normal_mode_distribution, normal_mode_transformation)
end

function Base.rand(rng::AbstractRNG, s::HarmonicRingPolymer)
    R = rand(rng, s.normal_mode_distribution)
    mul!(s.normal_mode_transformation.tmp, s.normal_mode_transformation.U, R)
    copy!(R, s.normal_mode_transformation.tmp)
    return R
end

function select_item(x::HarmonicRingPolymer{T}, ::Integer, size::NTuple) where {T}
    length(size) == 3 || throw(error("`size` $size must have 3 dimensions."))
    size[3] == length(x.normal_mode_distribution) || throw(error("Sample size does not match distribution."))

    out = Array{T,3}(undef, size)
    @views for i in axes(out, 2)
        for j in axes(out, 1)
            copyto!(out[j,i,:], rand(x))
        end
    end

    return out
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
