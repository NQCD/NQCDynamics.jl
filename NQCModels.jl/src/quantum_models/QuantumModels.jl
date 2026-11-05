
"""
    QuantumModels

Models defined within this module subtype the [`QuantumModel`](@ref) and provide
potentials as Hermitian matrices and derivatives as arrays of Hermitian matrices.
"""
module QuantumModels

using ..NQCModels
using NQCModels.BathDiscretisations: fillbathstates!, fillbathcoupling!, setcoupling!
using Unitful: @u_str
using UnitfulAtomic: austrip

using Parameters: Parameters
using LinearAlgebra: LinearAlgebra, Hermitian
using StaticArrays: SMatrix, SVector

"""
    QuantumModel <: Model

`QuantumModel`s are used when a system has multiple electronic states and the dynamics
of the system are propagated by a Hamiltonian in the diabatic representation. This is
the case for the majority of model systems.

# Implementation

`R` is always a `dimensions * atoms` shaped `AbstractMatrix`, `potential(model, R)` always
returns a single `nstates * nstates` Hermitian matrix (the sum of every atom's independent
contribution), and `derivative(model, R)`/`derivative!(model, D, R)` always return/fill a
`dimensions * atoms` shaped `Matrix` of `nstates * nstates` Hermitian matrices, one per
(degree of freedom, atom) pair.

Every model must implement:
* `potential!(model, V::Hermitian, R::AbstractMatrix)`: fill `V` with the potential for every
  atom (column) in `R`. Must fully overwrite `V.data` (e.g. zero it, then accumulate per atom) —
  never rely on `V` already being zero, since callers reuse the same array across calls.
* `derivative!(model, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)`: fill `D` in place
  with the derivative for every (degree of freedom, atom) pair in `R`.
* `nstates(model)`
* `ndofs(model)`

There is no automatic per-atom fallback: if every atom experiences the model
*independently* (the common case, with `ndofs(model) == 1`), the implementation should say so
explicitly, typically by looping over `R`/`axes(R, 2)` internally and accumulating each atom's
contribution. Models with `ndofs(model) > 1` (like [`GatesHollowayElbow`](@ref)), or that
explicitly implement an interaction between different atoms (like [`SpinBoson`](@ref)), are
written the same way — there is only one method to implement, not a separate "independent
atoms" special case.

The only automatic convenience is at the [`Model`](@ref NQCModels.Model) level: for models
where `ndofs(model) == 1`, the position of a single 1D atom may be passed directly as a plain
number rather than a `1×1` matrix.

# Example

In this example we create a simple 2 state, 1 dimensional quantum model `MyModel`, where
every atom contributes independently, then evaluate it for a single atom.

```jldoctest
using LinearAlgebra: Hermitian

struct MyModel <: NQCModels.QuantumModels.QuantumModel end

NQCModels.nstates(::MyModel) = 2
NQCModels.ndofs(::MyModel) = 1

function NQCModels.potential!(::MyModel, V::Hermitian, R::AbstractMatrix)
    V.data .= [sum(R) length(R); length(R) -sum(R)]
end

function NQCModels.derivative!(::MyModel, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    for atom in axes(R, 2)
        D[1,atom].data .= [1 0; 0 1]
    end
end

model = MyModel()
NQCModels.potential(model, hcat(10))

# output

2×2 Hermitian{Int64, Matrix{Int64}}:
 10    1
  1  -10
```
"""
abstract type QuantumModel <: NQCModels.Model end

"""
    QuantumFrictionModel <: QuantumModel

These models are defined identically to a typical `QuantumModel` but
allocate extra temporary arrays when used with `NQCDynamics.jl`.

This allows for the calculation of electronic friction
internally from the diabatic potential after diagonalisation
and calculation of nonadiabatic couplings.

When a molecular dynamics with electronic friction simulation is
set up in NQCDynamics, the `QuantumFrictionModel` is paired with a
`FrictionEvaluationMethod` in order to calculate the electronic
friction from the potential and derivative matrices.
"""
abstract type QuantumFrictionModel <: QuantumModel end

function matrix_template(model::QuantumModel, eltype)
    return zeros(eltype, NQCModels.nstates(model), NQCModels.nstates(model))
end

function vector_template(model::QuantumModel, eltype)
    return zeros(eltype, NQCModels.nstates(model))
end

function NQCModels.zero_derivative(model::QuantumModel, R::AbstractMatrix)
    return [Hermitian(matrix_template(model, eltype(R))) for _ in CartesianIndices(R)]
end

# GPU extensions: NQCModelsCUDAExt, NQCModelsMetalExt
# For GPU arrays, these return plain matrices rather than Hermitian-wrapped types
function NQCModels.potential(model::QuantumModel, R::AbstractMatrix)
    V = Hermitian(matrix_template(model, eltype(R)))
    NQCModels.potential!(model, V, R)
    return V
end

# GPU extensions: NQCModelsCUDAExt, NQCModelsMetalExt
# For GPU arrays, these return plain matrices rather than Hermitian-wrapped types
function NQCModels.derivative(model::QuantumModel, R::AbstractMatrix)
    D = NQCModels.zero_derivative(model, R)
    NQCModels.derivative!(model, D, R)
    return D
end

include("avoided_crossing_models/ananth_models.jl")
export AnanthModelOne
export AnanthModelTwo

include("avoided_crossing_models/tully_models.jl")
export TullyModelOne
export TullyModelTwo
export TullyModelThree

include("explicit_bath_models/anderson_holstein.jl")
export AndersonHolstein

include("explicit_bath_models/widebandbath.jl")
export WideBandBath

include("double_well.jl")
export DoubleWell

include("spin_boson.jl")
export OhmicSpectralDensity
export DebyeSpectralDensity
export SpinBoson
export BosonBath # Is a ClassicalModel type btw, just to really mess with the type structure.

include("erpenbeck_thoss.jl")
export ErpenbeckThoss

include("anderson_haldane.jl")
export AndersonHaldane

include("gates_holloway_elbow.jl")
export GatesHollowayElbow

include("subotnik.jl")
export MiaoSubotnik

include("three_state_morse.jl")
export ThreeStateMorse

include("1D_scattering.jl")
export Scattering1D

include("ouyang_models.jl")
export OuyangModelOne

include("adiabatic_state_selector.jl")
export AdiabaticStateSelector

end # module
