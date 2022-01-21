
module ClassicalMethods

using NonadiabaticMolecularDynamics:
    NonadiabaticMolecularDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsUtils,
    DynamicsMethods,
    Estimators,
    RingPolymers
using NonadiabaticMolecularDynamics.Calculators: Calculators
using NonadiabaticDynamicsBase: Atoms
using NonadiabaticModels: Model, ndofs

include("classical.jl")
export Classical
include("langevin.jl")
export Langevin, ThermalLangevin
include("mdef.jl")
export MDEF

const ClassicalMethodUnion = Union{Classical, Langevin, ThermalLangevin, MDEF}

function DynamicsUtils.classical_hamiltonian(sim::AbstractSimulation{<:ClassicalMethodUnion}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    potential = DynamicsUtils.classical_potential_energy(sim, DynamicsUtils.get_positions(u))
    return kinetic + potential
end

end # module
