module NQCDynamics

using Reexport: @reexport

@reexport using NQCBase
@reexport using NQCModels

include("RingPolymers/RingPolymers.jl")
@reexport using .RingPolymers
@reexport using RingPolymerArrays

include("Calculators/Calculators.jl")

include("simulations.jl")
export Simulation,
       RingPolymerSimulation,
       natoms,
       masses

include("NonadiabaticDistributions/NonadiabaticDistributions.jl")
@reexport using .NonadiabaticDistributions

include("DynamicsUtils/DynamicsUtils.jl")
@reexport using .DynamicsUtils: get_positions, get_velocities
export DynamicsUtils

include("Estimators.jl")
export Estimators

include("DynamicsOutputs.jl")

include("TimeCorrelationFunctions.jl")
export TimeCorrelationFunctions

include("DynamicsMethods/DynamicsMethods.jl")
@reexport using .DynamicsMethods

include("InitialConditions/InitialConditions.jl")
export InitialConditions

include("Ensembles/Ensembles.jl")
export Ensembles
@reexport using .Ensembles

end # module
