
module IntegrationAlgorithms

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsMethods,
    DynamicsUtils,
    natoms, nbeads, ndofs
using NQCDynamics.DynamicsMethods.SurfaceHoppingMethods
using NQCDynamics.DynamicsMethods.EhrenfestMethods
using NQCCalculators
using OrdinaryDiffEqCore: OrdinaryDiffEqCore, get_fsalfirstlast, OrdinaryDiffEqAlgorithm
using StochasticDiffEq: StochasticDiffEq

include("UniversalIntegrator.jl")
include("LegacyIntegrators/mdef_baoab.jl")
include("RingPolymerIntegrators/bcocb.jl")
include("mint.jl")
include("RingPolymerIntegrators/ringpolymer_mint.jl")
include("RingPolymerIntegrators/bcb.jl")
include("steps.jl")
include("LegacyIntegrators/verlet_with_electronics.jl")

#= include("CoupledIntegrators/CoupledIntegrator.jl")
include("CoupledIntegrators/create_problem.jl") =#

DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.ClassicalMethods.Classical}) = BAXAB(NoiseFree(), DynamicsMethods.ClassicalMethods.step_X_classical!)
DynamicsMethods.select_algorithm(sim::Simulation{<:DynamicsMethods.ClassicalMethods.LangevinMethods}) = BAXAB(NoiseCoupled(), DynamicsMethods.ClassicalMethods.step_X_constantfriction!, ConstantFriction(sim.method.γ))
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}, X::Function=noisedriven_step_X!) = BAXAB(NoiseCoupled(), X)#StochasticDiffEq.BAOAB(noise_mtx=true)#
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.MappingVariableMethods.SpinMappingW}) = MInt()
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping}, X::Function=electrondriven_step_X!) = BAXAB(NoiseFree(), X)
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.EhrenfestMethods.AbstractEhrenfest}, X::Function=electrondriven_step_X!) = BAXAB(NoiseFree(), X)

#RingPolymerSimulation Integrators

    # Classical Integrators
DynamicsMethods.select_algorithm(::RingPolymerSimulation{DynamicsMethods.ClassicalMethods.Classical}) = BCB()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.AbstractIESH}) = BCBWavefunction()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.ClassicalMasterEquation}) = BCBFull()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}) = BCOCB()

    # Ehrenfest methods
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.EhrenfestMethods.AbstractEhrenfest}) = BCBwithTsit5(OrdinaryDiffEq.Tsit5())
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.EhrenfestMethods.EhrenfestNA}) = BCBWavefunction()

    # MDEF methods
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}) = BCOCB()

    # Mapping Variable Integrators
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.MappingVariableMethods.eCMM}) = RingPolymerMInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.MappingVariableMethods.NRPMD}) = RingPolymerMInt()

    # Surface Hopping Methods
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping}) = BCBwithTsit5(OrdinaryDiffEq.Tsit5())

export BAXAB

end # module
