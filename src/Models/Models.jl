"""
    Models

Definitions for potentials and derivatives used to perform dynamics and sampling.

# Implementation
When implementing a new model, it is necessary to create a subtype of one of the child
types of `Model` and implement the functions required for the chosen type.
The concrete subtype can contain any parameters needed to evaluate the functions.
"""
module Models

using Unitful
using UnitfulAtomic
using LinearAlgebra
using Reexport
using Requires

using ..NonadiabaticMolecularDynamics

export Model
export AdiabaticModel
export DiabaticModel
export DiabaticFrictionModel
export FrictionModel

export potential!
export derivative!
export friction!

"""
    Model

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
    AdiabaticModel <: Model

This model implements `potential!` and `derivative!`.

`potential!` must fill an `AbstractVector` with `length = 1`.

`derivative!` must fill an `AbstractMatrix` with `size = (DoFs, atoms)`.
"""
abstract type AdiabaticModel <: Model end

"""
    DiabaticModel <: Model

This model implements `potential!` and `derivative!`.

`potential!` must fill a `Hermitian` with `size = (n_states, n_states)`.

`derivative!` must fill an `AbstractMatrix{<:Hermitian}` with `size = (DoFs, atoms)`,
and each entry must have `size = (n_states, n_states)`.
"""
abstract type DiabaticModel <: Model end

"""
    DiabaticFrictionModel <: Model

This model implements `potential!` and `derivative!`.

`potential!` must fill a `Hermitian` with `size = (n_states, n_states)`.

`derivative!` must fill an `AbstractMatrix{<:Hermitian}` with `size = (DoFs, atoms)`,
and each entry must have `size = (n_states, n_states)`.

Friction is calculated internally from the diabatic potential after diagonalisation
and calculation of nonadiabatic couplings.
"""
abstract type DiabaticFrictionModel <: DiabaticModel end

"""
    AdiabaticModel <: Model

This model implements `potential!`, `derivative!`, and `friction!`

`potential!` and `friction!` should be the same as for the `AdiabaticModel`.

`friction!` must fill an `AbstractMatrix` with `size = (DoFs*atoms, DoFs*atoms)`.
"""
abstract type FrictionModel <: Model end

"""
    potential!(model::Model, V, R::AbstractMatrix)

Fill `V` with the electronic potential of the system as a function of the positions `R`.

This must be implemented for all models.
"""
function potential! end

"""
    derivative!(model::Model, D, R::AbstractMatrix)

Fill `D` with the derivative of the electronic potential as a function of the positions `R`.

This must be implemented for all models.
"""
function derivative! end

"""
    friction!(model::FrictionModel, F, R:AbstractMatrix)

Fill `F` with the electronic friction as a function of the positions `R`.

This need only be implemented for `FrictionModel`s.
"""
function friction! end

include("analytic_models/free.jl")
include("analytic_models/harmonic.jl")
include("analytic_models/diatomic_harmonic.jl")

include("analytic_models/double_well.jl")
include("analytic_models/tully_models.jl")
include("analytic_models/scattering_anderson_holstein.jl")
include("analytic_models/1D_scattering.jl")

include("analytic_models/friction_harmonic.jl")
include("analytic_models/free_constant_friction.jl")

include("analytic_models/subotnik_A.jl")

include("EANN/EANN_H2Cu.jl")
include("EANN/EANN_H2Ag.jl")
include("EANN/EANN_NOAu.jl")

function __init__()
    @require PyCall="438e738f-606a-5dbb-bf0a-cddfbfd45ab0" @eval include("ML/ML.jl")
    @require JuLIP="945c410c-986d-556a-acb1-167a618e0462" @eval include("julip.jl")
end

function energy(model::Union{AdiabaticModel, FrictionModel}, R::AbstractMatrix)
    energy = [0.0]
    potential!(model, energy, R)
    energy[1]
end

function energy(model::DiabaticModel, R::AbstractMatrix; state=1)
    # The Diabatic Models contain a diabatic matrix
    V = Hermitian(zeros(model.n_states,model.n_states))
    # Get diabatic matrix
    potential!(model,V,R)
    # Export the eigenvalue lowest in energy
    eigvals(V)[state]
end

function eigen_summary(model::DiabaticModel, R::AbstractMatrix)
    # define output array
    eig_array = zeros(model.n_states+1, model.n_states)
    # The Diabatic Models contain a diabatic matrix
    V = Hermitian(zeros(model.n_states,model.n_states))
    # Get diabatic matrix
    potential!(model,V,R)
    # Export a matrix of eigenvectors and eigenvalues in the bottom line
    eig_array[1:model.n_states,:] .= eigvecs(V)
    eig_array[end,:] .= eigvals(V)
    eig_array = eig_array
end

function impurity_summary(model::DiabaticModel, R::AbstractMatrix, state::AbstractArray, σ::AbstractArray)
    # define output array
    #eig_array = zeros(model.n_states+2)
    eig_vec = zeros(model.n_states,model.n_states)
    eival = zeros(model.n_states)
    σdia = zeros(model.n_states,model.n_states)
    eig_array = zeros(4)
    # The Diabatic Models contain a diabatic matrix
    V = Hermitian(zeros(model.n_states,model.n_states))
    # Get diabatic matrix
    potential!(model,V,R)
    # Get the eigenvectors and values
    eival .= eigvals(V)
    eig_vec .= eigvecs(V)
    ieig = inv(eig_vec)
    #eig_array[1:model.n_states] .= eigvals(V)
    # calculate diabatic density matrix
    σdia .= eig_vec .* real(σ) .* ieig
    # Calculate two different hopping probabilities

    # save position
    eig_array[1] = R[1]
    for i = 1:length(state)
        # Energy
        eig_array[2] = eig_array[2] + state[i]*eival[i]
        # Hopping prob. by hopping array
        eig_array[3] = eig_array[3] + state[i]
        # directly from density matrix hopping prob.
        eig_array[4] = eig_array[4] + (real(σdia[i,i]) + imag(σdia[i,i]))^2*state[i]
    end
    #println(eig_array[2])
    # Export an array of eigenvalues with last two elements being hopping prob
    eig_array = eig_array
end

function potential_matrix(model::DiabaticModel, R::AbstractMatrix; state=1)
    # The Diabatic Models contain a diabatic matrix
    V = Hermitian(zeros(model.n_states,model.n_states))
    # Get diabatic matrix
    potential!(model,V,R)
    # Export the eigenvalue lowest in energy
    V = V
end

function forces(model::Union{AdiabaticModel, FrictionModel}, R::AbstractMatrix)
    forces = zero(R)
    derivative!(model, forces, R)
    -forces
end

include("plot.jl")
end # module
