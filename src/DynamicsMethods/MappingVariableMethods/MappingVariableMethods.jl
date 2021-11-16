
module MappingVariableMethods

using NonadiabaticMolecularDynamics:
    NonadiabaticMolecularDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsMethods,
    DynamicsUtils,
    Calculators,
    Estimators,
    NonadiabaticDistributions,
    TimeCorrelationFunctions,
    ndofs
using NonadiabaticModels: NonadiabaticModels, Model
using NonadiabaticDynamicsBase: Atoms
using .DynamicsUtils: get_mapping_momenta, get_mapping_positions

include("nrpmd.jl")
export NRPMD

include("cmm.jl")
export eCMM

end # module
