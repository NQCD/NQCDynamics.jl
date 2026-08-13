
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

# Implementation

`ClassicalModel`s should implement:
* `potential(model, R)`
* `derivative!(model, D, R)` (this is the derivative of the potential energy with respect to the positions)
* `ndofs(model)` (these are the degrees of freedom)

# Example

This example creates a 2 dimensional Classical model `MyModel`.
We implement the 3 compulsory functions then evaluate the potential.
Here, the argument `R` is an `AbstractMatrix` since this is a 2D model
that can accept multiple atoms.

```jldoctest
struct MyModel{P} <: NQCModels.ClassicalModels.ClassicalModel
    param::P
end

NQCModels.ndofs(::MyModel) = 2

NQCModels.potential(model::MyModel, R::AbstractMatrix) = model.param*sum(R.^2)
NQCModels.derivative!(model::MyModel, D, R::AbstractMatrix) = D .= model.param*2R

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
