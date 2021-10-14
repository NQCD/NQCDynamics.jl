
module RingPolymers

include("ring_polymer_array.jl")
export RingPolymerArray,
       get_centroid

include("ring_polymer.jl")
export transform_to_normal_modes!,
       transform_from_normal_modes!,
       nbeads

end # module
