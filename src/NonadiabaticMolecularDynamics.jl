module NonadiabaticMolecularDynamics

using Reexport
using Requires
using DocStringExtensions

include("unit_conversions.jl")

include("atoms.jl")
include("ring_polymer.jl")
include("ring_polymer_array.jl")
include("cells.jl")
include("dynamical_variables.jl")

include("Models/Models.jl")
export Models
include("Calculators/Calculators.jl")
export Calculators

include("simulations.jl")
include("classical_hamiltonians.jl")

include("InitialConditions/InitialConditions.jl")
export InitialConditions

include("Dynamics/Dynamics.jl")
@reexport using .Dynamics

function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" @eval include("ase_io.jl")
end

include("simulation_constructors.jl")

end # module
