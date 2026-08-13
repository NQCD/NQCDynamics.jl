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
"""
function potential(model::Model, R::AbstractMatrix)
    if ndofs(model) == 1
        if size(R, 2) == 1
            return potential(model, R[1])
        else
            return potential(model, @view R[1,:])
        end
    elseif size(R, 2) == 1
        return potential(model, @view R[:,1])
    else
        throw(MethodError(potential, (model, R)))
    end
end

"""
    potential(model::Model, R::Real)

Wraps R in a 1x1 matrix and redirects to potential(model::Model, R::AbstractMatrix).
"""
function potential(model::Model, R::Real)
    potential(model::Model, hcat(R))
end

"""
    potential!(model::Model, V, R::AbstractMatrix)

In-place version of `potential`, used to implement more efficient dynamics. 
"""
function potential!(model::Model, V, R::AbstractMatrix) end

"""
    potential!(model::Model, V, R::Real)

Wraps R in a 1x1 matrix and redirects to potential!(model::Model, V, R::AbstractMatrix). 
"""
function potential!(model::Model, V, R::Real) 
    potential!(model::Model, V, hcat(R))
end

"""
    derivative(model::Model, R)

Allocating version of `derivative!`, this definition should be suitable for almost all models.

Implement `zero_derivative` to allocate an appropriate array then implement `derivative!`
to fill the array.
"""
function derivative(model::Model, R)
    D = zero_derivative(model, R)
    derivative!(model, D, R)
    return D
end

"""
    derivative(model::Model, R::Real)

Wraps R in a 1x1 matrix and redirects to derivative(model::Model, R::AbstractMatrix). 
"""
function derivative(model::Model, R::Real)
    derivative(model::Model, hcat(R))
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
    derivative!(model::Model, D, R::Real)

Wraps R in a 1x1 matrix and redirects to derivative!(model::Model, D, R::AbstractMatrix). 
"""
function derivative!(model::Model, D, R::Real) 
    derivative!(model::Model, D, hcat(R))
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
