
using RecursiveArrayTools: ArrayPartition, NamedArrayPartition
using ComponentArrays: ComponentVector
using ..RingPolymers: RingPolymers

get_positions(u::ComponentVector) = u.r
get_velocities(u::ComponentVector) = u.v
get_flat_positions(u::ComponentVector) = @view get_positions(u)[:]
get_flat_velocities(u::ComponentVector) = @view get_velocities(u)[:]

# These must be kept since the symplectic integrators (VelocityVerlet) in DiffEq use ArrayPartitions.
get_positions(u::ArrayPartition) = u.x[2]
get_velocities(u::ArrayPartition) = u.x[1]

get_positions(u::NamedArrayPartition) = u.r
get_velocities(u::NamedArrayPartition) = u.v

function get_quantum_subsystem end

get_mapping_positions(u::ComponentVector) = u.qmap
get_mapping_momenta(u::ComponentVector) = u.pmap
get_mapping_positions(u::ComponentVector, i) = @view get_mapping_positions(u)[:,i]
get_mapping_momenta(u::ComponentVector, i) = @view get_mapping_momenta(u)[:,i]
