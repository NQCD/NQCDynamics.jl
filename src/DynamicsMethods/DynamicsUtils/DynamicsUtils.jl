
module DynamicsUtils

using ....NonadiabaticMolecularDynamics: AbstractSimulation, Simulation, RingPolymerSimulation
using ..Dynamics: Dynamics

divide_by_mass!(dv, masses) = dv ./= masses'
velocity!(dr, v, r, sim, t) = dr .= v

"""
    apply_interbead_coupling!(du::DynamicalVariables, u::DynamicalVariables,
                              sim::RingPolymerSimulation)
    
Applies the force that arises from the harmonic springs between adjacent beads.

Only applies the force for atoms labelled as quantum within the `RingPolymerParameters`.
"""
function apply_interbead_coupling!(dr::AbstractArray{T,3}, r::AbstractArray{T,3}, sim::RingPolymerSimulation) where {T}
    for i in axes(dr, 3)
        iplus = mod1(i+1, sim.beads.n_beads)
        iminus = mod1(i-1, sim.beads.n_beads)
        for j in sim.beads.quantum_atoms
            for k in axes(dr, 1)
                dr[k,j,i] -= sim.beads.ω_n² * (2r[k,j,i] - r[k,j,iplus] - r[k,j,iminus])
            end
        end
    end
    return nothing
end

include("callbacks.jl")
export CellBoundaryCallback
export TerminatingCallback

include("output.jl")
include("density_matrix_dynamics.jl")
include("plot.jl")

end # module
