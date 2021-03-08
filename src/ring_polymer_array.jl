using RecursiveArrayTools
using StatsBase: mean
using LinearAlgebra: mul!

export RingPolymerArray
export transform!
export get_centroid

mutable struct RingPolymerArray{T} <: AbstractVectorOfArray{T,3,Vector{Matrix{T}}}
    u::Vector{Matrix{T}}
    U::Matrix{T}
    normal::Bool
    tmp::Vector{T}
    function RingPolymerArray{T}(u) where {T}
        U = get_normal_mode_transformation(length(u))
        new(u, U, false, zeros(length(u)))
    end
end

function RingPolymerArray(vec::Vector{T}) where {T}
    RingPolymerArray{eltype(T)}(vec)
end

"""
Creates normal mode transformation for `n` beads.
"""
function get_normal_mode_transformation(n::Integer)::Matrix

    # Real and imaginary parts of the discrete Fourier transform matrix. 
    U = sqrt(2/n) .* hcat([cos(2π * j * k / n) for j=0:n-1, k=0:n÷2],
                        [sin(2π * j * k / n) for j=0:n-1, k=n÷2+1:n-1])

    # Normalisation
    U[:, 1] ./= sqrt(2)
    iseven(n) && (U[:, n÷2+1] ./= sqrt(2))

    U
end

"""
Transform to/from normal mode coordinates.
"""
function transform!(z::RingPolymerArray)
    transform = z.normal ? z.U : z.U'
    @views for I in CartesianIndices(z[1])
        mul!(z.tmp, transform, z[I,:])
        z[I,:] .= z.tmp
    end
    z.normal = !z.normal
end

"""
Evaluate centroid of ring polymer.
"""
function get_centroid(z::RingPolymerArray{T})::Matrix{T} where {T}
    z.normal ? z[1] ./ sqrt(length(z)) : mean(z.u)
end
