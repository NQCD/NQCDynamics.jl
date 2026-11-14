"""
NQCModels define the potentials and derivatives that
govern the dynamics of the particles.
These can exist as analytic models or as interfaces to other codes.
"""
module NQCModels

using Reexport: @reexport

export potential, potential!
export derivative, derivative!
export nstates
export ndofs

"""
Top-level type for models.

# Implementation
When adding new models, this should not be directly subtyped. Instead, depending on
the intended functionality of the model, one of the child abstract types should be
subtyped.
If an appropriate type is not already available, a new abstract subtype should be created.
"""
abstract type Model end

Base.broadcastable(model::Model) = Ref(model)

"""
    potential(model::Model, R::AbstractMatrix)

Evaluate the potential at position `R` for the given `model`.

# GPU Support
GPU array support is provided via package extensions:
- `NQCModelsCUDAExt`: CUDA.jl support (activated by `using CUDA`)
- `NQCModelsMetalExt`: Metal.jl support (activated by `using Metal`)

When GPU arrays are passed, the return type matches the input array type.
"""
function potential(model::Model, R::AbstractMatrix)
    throw(MethodError(potential, (model, R)))
end

"""
    potential(model::Model, R::Number)

Convenience method for models where `ndofs(model) == 1`, allowing the position of a
single 1D atom to be passed directly as a number rather than a `1×1` matrix.
Wraps `R` in a `1×1` matrix and redirects to `potential(model::Model, R::AbstractMatrix)`.

Typed as `Number`, rather than `Real`, to stay unambiguous should a model ever dispatch
on a `Number`-typed argument itself.
"""
function potential(model::Model, R::Number)
    potential(model, hcat(R))
end


"""
    potential!(model::Model, V, R::AbstractMatrix)

In-place version of `potential`, used to implement more efficient dynamics.
"""
function potential!(model::Model, V, R::AbstractMatrix) end

"""
    potential!(model::Model, V, R::Number)

Convenience method for models where `ndofs(model) == 1`, allowing the position of a
single 1D atom to be passed directly as a number rather than a `1×1` matrix.
Wraps `R` in a `1×1` matrix and redirects to `potential!(model::Model, V, R::AbstractMatrix)`.

Typed as `Number`, rather than `Real`, to stay unambiguous should a model ever dispatch
on a `Number`-typed argument itself.
"""
function potential!(model::Model, V, R::Number)
    potential!(model, V, hcat(R))
end

"""
    derivative(model::Model, R)

Allocating version of `derivative!`, this definition should be suitable for almost all models.

Implement `zero_derivative` to allocate an appropriate array then implement `derivative!`
to fill the array.

# GPU Support
GPU array support is provided via package extensions:
- `NQCModelsCUDAExt`: CUDA.jl support (activated by `using CUDA`)
- `NQCModelsMetalExt`: Metal.jl support (activated by `using Metal`)

When GPU arrays are passed, the return type matches the input array type.
"""
function derivative(model::Model, R)
    D = zero_derivative(model, R)
    derivative!(model, D, R)
    return D
end

"""
    derivative(model::Model, R::Number)

Convenience method for models where `ndofs(model) == 1`, allowing the position of a
single 1D atom to be passed directly as a number rather than a `1×1` matrix.
Wraps `R` in a `1×1` matrix and redirects to `derivative(model::Model, R::AbstractMatrix)`.

Typed as `Number`, rather than `Real`, to stay unambiguous should a model ever dispatch
on a `Number`-typed argument itself.
"""
function derivative(model::Model, R::Number)
    derivative(model, hcat(R))
end


"""
    derivative!(model::Model, D, R::AbstractMatrix)

Fill `D` with the derivative of the electronic potential as a function of the positions `R`.

This must be implemented for all models.
"""
function derivative!(model::Model, D, R::AbstractMatrix)
    @warn "You have likely made a multiple dispatch mistake." maxlog = 1 model = model D = D R = R
end

"""
    derivative!(model::Model, D, R::Number)

Convenience method for models where `ndofs(model) == 1`, allowing the position of a
single 1D atom to be passed directly as a number rather than a `1×1` matrix.
Wraps `R` in a `1×1` matrix and redirects to `derivative!(model::Model, D, R::AbstractMatrix)`.

Typed as `Number`, rather than `Real`, to stay unambiguous should a model ever dispatch
on a `Number`-typed argument itself.
"""
function derivative!(model::Model, D, R::Number)
    derivative!(model, D, hcat(R))
end


"""
    zero_derivative(model::Model, R)

Create an zeroed array of the right size to match the derivative.
"""
function zero_derivative(model::Model, R) end

"""
    nstates(::Model)

Get the number of electronic states in the model.
"""
nstates(::Model) = error("This should return the total number of electronic states.")

"""
    ndofs(::Model)

Get the number of degrees of freedom for every atom in the model. Usually 1 or 3.
"""
ndofs(::Model) = error("This should return the number of degrees of freedom for each atom.")

dofs(model::Model) = Base.OneTo(ndofs(model))

function get_subsystem_derivative() end

state_independent_potential(model, r) = 0.0

function state_independent_potential!(model, Vsystem, r)
    Vsystem .= 0.0
end

function state_independent_derivative(model, r) end
state_independent_derivative!(model, derivative, r) = fill!(derivative, zero(eltype(r)))
nelectrons(::Model) = error("This should return the total number of electrons.")
fermilevel(::Model) = 0.0

eachelectron(model::Model) = Base.OneTo(nelectrons(model))
eachstate(model::Model) = Base.OneTo(nstates(model))

mobileatoms(::Model, n::Int) = Base.OneTo(n)

include("bath_discretisations/BathDiscretisations.jl")
@reexport using .BathDiscretisations

include("classical_models/ClassicalModels.jl")
@reexport using .ClassicalModels

include("quantum_models/QuantumModels.jl")
@reexport using .QuantumModels

include("friction_models/FrictionModels.jl")
@reexport using .FrictionModels

include("plot.jl")

include("subsystems.jl")

end # module
