module Models

export Model

abstract type Model end

include("Analytic/Analytic.jl")
include("ML/ML.jl")
end # module
