module NQCDynamics

using Reexport: @reexport

@reexport using NQCBase
@reexport using NQCModels

include("NumericUtils/FastDeterminant.jl")

include("RingPolymers/RingPolymers.jl")
@reexport using .RingPolymers
@reexport using RingPolymerArrays

include("Calculators/Calculators.jl")

include("simulations.jl")
export Simulation,
       RingPolymerSimulation,
       natoms,
       masses,
       Thermostat,
       get_temperature

@reexport using NQCDistributions
# Simulation-aware version of nuclear Boltzmann distribution. 
include("NQCDistributions-convenience.jl")

include("DynamicsUtils/DynamicsUtils.jl")
@reexport using .DynamicsUtils: get_positions, get_velocities
export DynamicsUtils

# Needs Simulation, NQCBase
include("structure.jl")
export Structure

include("Estimators.jl")
export Estimators

include("TimeCorrelationFunctions.jl")
export TimeCorrelationFunctions

include("DynamicsMethods/DynamicsMethods.jl")
@reexport using .DynamicsMethods

include("InitialConditions/InitialConditions.jl")
export InitialConditions

include("Ensembles/Ensembles.jl")
export Ensembles
@reexport using .Ensembles

# Needs DynamicsUtils, Simulation, InitialConditions
include("Analysis/Analysis.jl")
export Analysis

# Needs Analysis
include("DynamicsOutputs.jl")
@reexport using .DynamicsOutputs

end # module
