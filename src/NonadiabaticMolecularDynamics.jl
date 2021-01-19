module NonadiabaticMolecularDynamics

using Reexport
using Requires

include("atoms.jl")
include("ring_polymer.jl")
include("cells.jl")
include("phasespace.jl")

include("Models/Models.jl")
export Models
include("Calculators/Calculators.jl")
export Calculators

include("simulations.jl")
include("classical_hamiltonians.jl")

include("InitialConditions/InitialConditions.jl")
export InitialConditions

include("Dynamics/Dynamics.jl")
export Dynamics
using .Dynamics: SurfaceHoppingPhasespace
export SurfaceHoppingPhasespace

function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" @eval include("ase_io.jl")
end

end # module
