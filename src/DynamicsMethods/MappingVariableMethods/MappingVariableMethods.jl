
module MappingVariableMethods

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsMethods,
    DynamicsUtils,
    Estimators,
    TimeCorrelationFunctions,
    ndofs
using NQCCalculators
using NQCModels: NQCModels, Model
using NQCBase: Atoms
using .DynamicsUtils: get_mapping_momenta, get_mapping_positions

include("nrpmd.jl")
export NRPMD

include("cmm.jl")
export eCMM
include("rpecmm.jl")

include("spin_mapping.jl")
export SpinMappingW

end # module
