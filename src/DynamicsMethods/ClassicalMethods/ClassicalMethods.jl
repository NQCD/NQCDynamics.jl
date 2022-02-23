
module ClassicalMethods

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsUtils,
    DynamicsMethods,
    Estimators,
    RingPolymers
using NQCDynamics.Calculators: Calculators
using NQCBase: Atoms
using NQCModels: Model, ndofs

include("classical.jl")
export Classical
include("langevin.jl")
export Langevin, ThermalLangevin
include("mdef.jl")
export MDEF
include("diabatic_mdef.jl")
export DiabaticMDEF

const ClassicalMethodUnion = Union{Classical, Langevin, ThermalLangevin, MDEF, DiabaticMDEF}

function DynamicsUtils.classical_hamiltonian(sim::AbstractSimulation{<:ClassicalMethodUnion}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    potential = DynamicsUtils.classical_potential_energy(sim, DynamicsUtils.get_positions(u))
    return kinetic + potential
end

end # module
