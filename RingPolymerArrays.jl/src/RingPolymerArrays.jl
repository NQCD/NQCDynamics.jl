module RingPolymerArrays

using LinearAlgebra: normalize!, mul!
using OrderedCollections: OrderedDict

include("ring_polymer_array.jl")
export RingPolymerArray
export eachbead
export eachdof

include("normal_mode_transformation.jl")
export NormalModeTransformation
export transform_to_normal_modes!
export transform_from_normal_modes!

end # module
