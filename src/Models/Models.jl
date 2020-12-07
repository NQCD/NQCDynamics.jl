module Models

using Distances
using ..Atoms

export Model
export get_distances

abstract type Model end

function get_distances(positions::Vector, n_DoF::Integer=3)
    R = reshape(positions, n_DoF, :)
    pairwise(Euclidean(), R, dims=2)
end

include("Analytic/Analytic.jl")
include("ML/ML.jl")
include("EAM/PdH.jl")
end # module
