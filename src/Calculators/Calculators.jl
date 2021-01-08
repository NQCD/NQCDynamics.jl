module Calculators

using LinearAlgebra
using ..Models

abstract type AbstractCalculator{M<:Model} end
abstract type AbstractAdiabaticCalculator{M<:AdiabaticModel} <: AbstractCalculator{M} end
abstract type AbstractDiabaticCalculator{M<:DiabaticModel} <: AbstractCalculator{M} end
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
    function DiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        potential = Hermitian(zeros(model.n_states, model.n_states))
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms]
        eigenvalues = zeros(model.n_states)
        eigenvectors = zeros(model.n_states, model.n_states)
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors, adiabatic_derivative)
    end
end

struct RingPolymerDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Vector{Hermitian{T}}
    derivative::Array{Hermitian{T},3}
    eigenvalues::Vector{Vector{T}}
    eigenvectors::Vector{Matrix{T}}
    adiabatic_derivative::Array{Matrix{T},3}
    function RingPolymerDiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        potential = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:beads]
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms, k=1:beads]
        eigenvalues = [zeros(model.n_states) for i=1:beads]
        eigenvectors = [zeros(model.n_states, model.n_states) for i=1:beads]
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms, k=1:beads]
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors, adiabatic_derivative)
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

function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    DiabaticCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    AdiabaticCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::FrictionModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    FrictionCalculator{T}(model, DoFs, atoms)
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

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractMatrix{T}) where {T}
    Models.friction!(calc.model, calc.friction, R)
end

function eigen!(calc::DiabaticCalculator)
    eig = eigen(calc.potential)
    @views for i=1:length(eig.values)
        if eig.vectors[:,i]'calc.eigenvectors[:,i] < 0
            calc.eigenvectors[:,i] .= -eig.vectors[:,i]
        else
            calc.eigenvectors[:,i] .= eig.vectors[:,i]
        end
    end
    calc.eigenvalues .= eig.values
end

function eigen!(calc::RingPolymerDiabaticCalculator)
    eigs = eigen.(calc.potential)
    calc.eigenvalues .= [eig.values for eig in eigs]
    calc.eigenvectors .= [eig.vectors for eig in eigs]
end

function transform_derivative!(calc::DiabaticCalculator)
    for i in axes(calc.derivative, 2) # Atoms
        for j in axes(calc.derivative, 1) # DoFs
            calc.adiabatic_derivative[j,i] .= calc.eigenvectors' * calc.derivative[j,i] * calc.eigenvectors
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

end # module