"""
    Calculators
This module exists to bridge the gap between the `Models` and the `Dynamics`.
Here we provide functions and types for evaluating and storing quantities obtained from the
`Models`.
In addition any further manipulation of those quantities, such as computing eigenvalues,
is included here.
This module is largely needed to facilitate seemless integration of both ring polymer and
classical dynamics, using the same models and functions for both.
Specific ring polymer types are provided that have the extra fields and methods needed
to evaluate the quantities for each bead. 
"""
module Calculators

using LinearAlgebra
using ..Models

"""
    AbstractCalculator{M<:Model}
Top-level type for all calculators.
Each concrete calculator contains the `Model` and the fields to store the quantities
obtained from the model.
"""
abstract type AbstractCalculator{M<:Model} end
abstract type AbstractAdiabaticCalculator{M<:AdiabaticModel} <: AbstractCalculator{M} end
abstract type AbstractDiabaticCalculator{M<:Union{DiabaticFrictionModel,DiabaticModel}} <: AbstractCalculator{M} end
abstract type AbstractFrictionCalculator{M<:FrictionModel} <: AbstractCalculator{M} end

struct AdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{M}
    model::M
    potential::Vector{T}
    derivative::Matrix{T}
    function AdiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, [0.0], zeros(DoFs, atoms))
    end
end

struct RingPolymerAdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{M}
    model::M
    potential::Vector{Vector{T}}
    derivative::Array{T,3}
    function RingPolymerAdiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        new{T,M}(model, [[0.0] for i=1:beads], zeros(DoFs, atoms, beads))
    end
end

struct DiabaticCalculator{T,M} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Hermitian{T}
    derivative::Matrix{Hermitian{T}}
    eigenvalues::Vector{T}
    eigenvectors::Matrix{T}
    adiabatic_derivative::Matrix{Matrix{T}}
    nonadiabatic_coupling::Matrix{Matrix{T}}
    tmp_mat1::Matrix{Complex{T}}
    tmp_mat2::Matrix{Complex{T}}
    tmp_mat3::Matrix{Complex{T}}
    function DiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        potential = Hermitian(zeros(model.n_states, model.n_states))
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms]
        eigenvalues = zeros(model.n_states)
        eigenvectors = zeros(model.n_states, model.n_states)
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        nonadiabatic_coupling = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        tmp_mat1 = zeros(model.n_states, model.n_states)
        tmp_mat2 = zeros(model.n_states, model.n_states)
        tmp_mat3 = zeros(model.n_states, model.n_states)
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors, 
                adiabatic_derivative, nonadiabatic_coupling, tmp_mat1,
                tmp_mat2, tmp_mat3)
    end
end

struct RingPolymerDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Vector{Hermitian{T}}
    derivative::Array{Hermitian{T},3}
    eigenvalues::Vector{Vector{T}}
    eigenvectors::Vector{Matrix{T}}
    adiabatic_derivative::Array{Matrix{T},3}
    nonadiabatic_coupling::Array{Matrix{T},3}
    function RingPolymerDiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        potential = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:beads]
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms, k=1:beads]
        eigenvalues = [zeros(model.n_states) for i=1:beads]
        eigenvectors = [zeros(model.n_states, model.n_states) for i=1:beads]
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms, k=1:beads]
        nonadiabatic_coupling = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms, k=1:beads]
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors, adiabatic_derivative, nonadiabatic_coupling)
    end
end

struct FrictionCalculator{T,M} <: AbstractFrictionCalculator{M}
    model::M
    potential::Vector{T}
    derivative::Matrix{T}
    friction::Matrix{T}
    function FrictionCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, [0.0], zeros(DoFs, atoms), zeros(DoFs*atoms, DoFs*atoms))
    end
end

struct DiabaticFrictionCalculator{T,M} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Hermitian{T}
    derivative::Matrix{Hermitian{T}}
    eigenvalues::Vector{T}
    eigenvectors::Matrix{T}
    adiabatic_derivative::Matrix{Matrix{T}}
    nonadiabatic_coupling::Matrix{Matrix{T}}
    friction::Matrix{T}
    function DiabaticFrictionCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        potential = Hermitian(zeros(model.n_states, model.n_states))
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms]
        eigenvalues = zeros(model.n_states)
        eigenvectors = zeros(model.n_states, model.n_states)
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        nonadiabatic_coupling = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        friction = zeros(DoFs*atoms, DoFs*atoms)
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors,
                 adiabatic_derivative, nonadiabatic_coupling, friction)
    end
end

function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    DiabaticCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    AdiabaticCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::FrictionModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    FrictionCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::DiabaticFrictionModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    DiabaticFrictionCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerDiabaticCalculator{T}(model, DoFs, atoms, beads)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerAdiabaticCalculator{T}(model, DoFs, atoms, beads)
end

function evaluate_potential!(calc::AbstractCalculator, R::AbstractMatrix{T}) where {T}
    Models.potential!(calc.model, calc.potential, R)
end

function evaluate_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        Models.potential!(calc.model, calc.potential[i], R[:,:,i])
    end
end

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractMatrix{T}) where {T}
    Models.derivative!(calc.model, calc.derivative, R)
end

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        Models.derivative!(calc.model, calc.derivative[:,:,i], R[:,:,i])
    end
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractMatrix)
    Models.friction!(calc.model, calc.friction, R)
end

@doc raw"""
    evaluate_friction!(calc::DiabaticFrictionCalculator, R::AbstractMatrix)
Evaluate the electronic friction for a model given in the diabatic representation.
Requires that `adiabatic_derivative` and `eigenvalues` be precomputed.
```math
γ = 2πħ ∑ⱼ <1|dH|j><j|dH|1> δ(ωⱼ) / ωⱼ
```
Note that the delta function is approximated by a normalised gaussian.
"""
function evaluate_friction!(calc::DiabaticFrictionCalculator, R::AbstractMatrix)

    gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))

    DoFs = size(R)[1]
    calc.friction .= 0
    for i in axes(R, 2) # Atoms
        for j in axes(R, 1) # DoFs
            for m=2:calc.model.n_states
                ω = calc.eigenvalues[m] - calc.eigenvalues[1]
                g = gauss(ω, calc.model.σ) / ω
                calc.friction[j+(i-1)*DoFs] += 2π*abs2(calc.adiabatic_derivative[j,i][m,1])*g
            end
        end
    end
end

function eigen!(calc::AbstractDiabaticCalculator)
    eig = eigen(calc.potential)
    correct_phase!(eig, calc.eigenvectors)
    copyto!(calc.eigenvectors, eig.vectors)
    copyto!(calc.eigenvalues, eig.values)
end

function eigen!(calc::RingPolymerDiabaticCalculator)
    eigs = eigen.(calc.potential)
    correct_phase!.(eigs, calc.eigenvectors)
    calc.eigenvalues .= [eig.values for eig in eigs]
    calc.eigenvectors .= [eig.vectors for eig in eigs]
end

function correct_phase!(eig::Eigen, old_eigenvectors::Matrix)
    @views for i=1:length(eig.values)
        if dot(eig.vectors[:,i], old_eigenvectors[:,i]) < 0
            eig.vectors[:,i] .*= -1
        end
    end
end

# function transform_derivative!(calc::AbstractDiabaticCalculator)
#     for i in axes(calc.derivative, 2) # Atoms
#         for j in axes(calc.derivative, 1) # DoFs
#             calc.adiabatic_derivative[j,i] .= calc.eigenvectors' * calc.derivative[j,i] * calc.eigenvectors
#         end
#     end
# end

function transform_derivative!(calc::AbstractDiabaticCalculator)
    # speed up calculation with mul()
    calc.tmp_mat2 .= calc.eigenvectors
    for i in axes(calc.derivative, 2) # Atoms
        for j in axes(calc.derivative, 1) # DoFs
            calc.tmp_mat3 .= calc.derivative[j,i]
            mul!(calc.tmp_mat1, calc.tmp_mat3, calc.tmp_mat2)
            mul!(calc.tmp_mat3, calc.tmp_mat2', calc.tmp_mat1)
            calc.adiabatic_derivative[j,i] .= calc.tmp_mat3
        end
    end
end

function transform_derivative!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.derivative, 3) # Beads
        for j in axes(calc.derivative, 2) # Atoms
            for k in axes(calc.derivative, 1) # DoFs
                calc.adiabatic_derivative[k,j,i] .= calc.eigenvectors[i]' * calc.derivative[k,j,i] * calc.eigenvectors[i]
            end
        end
    end
end


function evaluate_nonadiabatic_coupling!(calc::DiabaticCalculator)
    evaluate_nonadiabatic_coupling!.(calc.nonadiabatic_coupling, calc.adiabatic_derivative,
                                     Ref(calc.eigenvalues))
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.nonadiabatic_coupling, 3) # Beads
        evaluate_nonadiabatic_coupling!.(calc.nonadiabatic_coupling[:,:,i], calc.adiabatic_derivative[:,:,i],
                                        Ref(calc.eigenvalues[i]))
    end
end

function evaluate_nonadiabatic_coupling!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalues::Vector)
    for i=1:length(eigenvalues)
        for j=i+1:length(eigenvalues)
            coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
            coupling[i,j] = -coupling[j,i]
        end
    end
end

function update_electronics!(calculator::AbstractDiabaticCalculator, r::AbstractMatrix)
    evaluate_potential!(calculator, r)
    evaluate_derivative!(calculator, r)
    eigen!(calculator)
    transform_derivative!(calculator)
    evaluate_nonadiabatic_coupling!(calculator)
end

end # module