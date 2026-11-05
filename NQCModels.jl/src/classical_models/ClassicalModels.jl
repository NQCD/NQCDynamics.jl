
"""
    ClassicalModels

All models defined within this module have only a single electronic state
and return potentials as scalars and derivatives as simple arrays.

The central abstract type is the [`ClassicalModel`](@ref), which all models should subtype.
"""
module ClassicalModels

using ..NQCModels: NQCModels
using Requires: Requires
using Parameters: Parameters

using Unitful: @u_str, ustrip
using UnitfulAtomic: austrip, auconvert

using Reexport: @reexport

"""
    ClassicalModel <: Model

`ClassicalModel`s represent the potentials from classical molecular dynamics
where the potential is a function of the position.

`R` is always a `dimensions * atoms` shaped `AbstractMatrix`, `potential(model, R)`
always returns a scalar, and `derivative(model, R)`/`derivative!(model, D, R)` always
return/fill a `dimensions * atoms` shaped array of scalars.

# Implementation

Every model must implement:
* `potential(model, R::AbstractMatrix)`: the total potential for every atom (column) in `R`.
* `derivative!(model, D::AbstractMatrix, R::AbstractMatrix)`: fill `D` in place with the
  derivative for every atom (column) in `R`.
* `ndofs(model)` (the degrees of freedom per atom)

There is no automatic per-atom fallback: if every atom experiences the model
*independently*, the implementation should say so explicitly, typically by summing/broadcasting
over `eachcol(R)`/`axes(R, 2)` internally. Models that explicitly implement an interaction
between different atoms (like [`DiatomicHarmonic`](@ref)), or that support heterogeneous
per-atom parameters (like [`Harmonic`](@ref)), are written the same way — there is only one
method to implement, not a separate "independent atoms" special case.

The only automatic convenience is at the [`Model`](@ref NQCModels.Model) level: for models
where `ndofs(model) == 1`, the position of a single 1D atom may be passed directly as a plain
number rather than a `1×1` matrix.

# Example

This example creates a 2 dimensional Classical model `MyModel` where every atom is
independent, then evaluates the potential for a 2 atom system.

```jldoctest
struct MyModel{P} <: NQCModels.ClassicalModels.ClassicalModel
    param::P
end

NQCModels.ndofs(::MyModel) = 2

NQCModels.potential(model::MyModel, R::AbstractMatrix) = model.param*sum(R.^2)
function NQCModels.derivative!(model::MyModel, D::AbstractMatrix, R::AbstractMatrix)
    D .= model.param*2R
end

model = MyModel(10)

NQCModels.potential(model, [1 2; 3 4])

# output

300

```
"""
abstract type ClassicalModel <: NQCModels.Model end

NQCModels.nstates(::ClassicalModel) = 1

NQCModels.zero_derivative(::ClassicalModel, R) = zero(R)

function NQCModels.potential!(model::ClassicalModel, V::Real, R::AbstractMatrix)
    @warn "In Julia, real numbers cannot be updated in place. Please define your potential as a matrix if you wish to pass it into potential!()"
end

"""
    NQCModels.potential!(model::ClassicalModel, V::AbstractMatrix, R::AbstractMatrix)

Default in-place `potential!`, suitable for every `ClassicalModel` since `potential`
always returns a scalar.
"""
function NQCModels.potential!(model::ClassicalModel, V::AbstractMatrix, R::AbstractMatrix)
    V .= NQCModels.potential(model, R)
end

# GPU extensions: NQCModelsCUDAExt, NQCModelsMetalExt
function NQCModels.derivative(model::ClassicalModel, R::AbstractMatrix)
    D = NQCModels.zero_derivative(model, R)
    NQCModels.derivative!(model, D, R)
    return D
end

include("free.jl")
export Free
include("morse.jl")
export Morse
include("logistic.jl")
export Logistic
include("harmonic.jl")
export Harmonic
include("diatomic_harmonic.jl")
export DiatomicHarmonic
include("darling_holloway_elbow.jl")
export DarlingHollowayElbow
include("averaged_potential.jl")
export AveragedPotential
include("csv_models.jl")
export CSVModel_1D
include("atomscalculators.jl")
export AtomsCalculatorsModel

function __init__()
    Requires.@require JuLIP="945c410c-986d-556a-acb1-167a618e0462" @eval include("julip.jl")
end

end # module
