
module IntegrationAlgorithms

using NQCDynamics:
    NQCDynamics,
    AbstractSimulation,
    Simulation,
    RingPolymerSimulation,
    DynamicsMethods,
    DynamicsUtils,
    Calculators,
    natoms, nbeads, ndofs

using OrdinaryDiffEq: OrdinaryDiffEq
using StochasticDiffEq: StochasticDiffEq

struct BCB <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end
struct BCBwithTsit5 <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

"""
    RingPolymerMInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm

Second order symplectic momentum integral algorithm applied to NRPMD.

# Reference

[J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
"""
struct RingPolymerMInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

"""
    MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm

Second order symplectic momentum integral algorithm.

# Reference

[J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
"""
struct MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

struct MDEF_BAOAB <: StochasticDiffEq.StochasticDiffEqAlgorithm end
struct BCOCB <: StochasticDiffEq.StochasticDiffEqAlgorithm end

DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.SurfaceHoppingMethods.SurfaceHopping}) = BCBwithTsit5()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.EhrenfestMethods.AbstractEhrenfest}) = BCBwithTsit5()
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}) = MDEF_BAOAB()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}) = BCOCB()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}) = BCOCB()
DynamicsMethods.select_algorithm(::Simulation{<:DynamicsMethods.MappingVariableMethods.SpinMappingW}) = MInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.MappingVariableMethods.NRPMD}) = RingPolymerMInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.MappingVariableMethods.eCMM}) = RingPolymerMInt()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{DynamicsMethods.ClassicalMethods.Classical}) = BCB()

export BCB
export BCBwithTsit5
export RingPolymerMInt
export MInt
export MDEF_BAOAB
export BCOCB

include("mdef_baoab.jl")
include("bcocb.jl")
include("mint.jl")
include("ringpolymer_mint.jl")
include("bcb_electronics.jl")
include("bcb.jl")
include("steps.jl")

end # module
