
using RecursiveArrayTools: ArrayPartition
using ComponentArrays: ComponentVector

get_positions(u::ComponentVector) = u.r
get_velocities(u::ComponentVector) = u.v
get_flat_positions(u::ComponentVector) = @view get_positions(u)[:]
get_flat_velocities(u::ComponentVector) = @view get_velocities(u)[:]

# These must be kept since the symplectic integrators (VelocityVerlet) in DiffEq use ArrayPartitions.
get_positions(u::ArrayPartition) = u.x[2]
get_velocities(u::ArrayPartition) = u.x[1]

function get_quantum_subsystem end
