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
using NonadiabaticModels
using ..NonadiabaticMolecularDynamics: get_centroid
using StaticArrays

"""
    AbstractCalculator{M<:Model}

Top-level type for all calculators.

Each concrete calculator contains the `Model` and the fields to store the quantities
obtained from the model.
"""
abstract type AbstractCalculator{M<:Model} end
abstract type AbstractAdiabaticCalculator{M<:AdiabaticModel} <: AbstractCalculator{M} end
abstract type AbstractDiabaticCalculator{M<:Union{DiabaticFrictionModel,DiabaticModel}} <: AbstractCalculator{M} end
abstract type AbstractStaticDiabaticCalculator{M} <: AbstractDiabaticCalculator{M} end
abstract type AbstractFrictionCalculator{M<:AdiabaticFrictionModel} <: AbstractCalculator{M} end

mutable struct AdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{M}
    model::M
    potential::T
    derivative::Matrix{T}
    function AdiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, 0, zeros(DoFs, atoms))
    end
end

struct RingPolymerAdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{M}
    model::M
    potential::Vector{T}
    derivative::Array{T,3}
    function RingPolymerAdiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        new{T,M}(model, zeros(beads), zeros(DoFs, atoms, beads))
    end
end

mutable struct DiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Hermitian{T,SMatrix{S,S,T,L}}
    derivative::Matrix{Hermitian{T,SMatrix{S,S,T,L}}}
    eigenvalues::SVector{S,T}
    eigenvectors::SMatrix{S,S,T,L}
    adiabatic_derivative::Matrix{SMatrix{S,S,T,L}}
    nonadiabatic_coupling::Matrix{SMatrix{S,S,T,L}}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function DiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        n = Int(model.n_states)
        matrix_template = NonadiabaticModels.matrix_template(model, T)
        vector_template = NonadiabaticModels.vector_template(model, T)

        potential = Hermitian(matrix_template)
        derivative = [Hermitian(matrix_template) for _=1:DoFs, _=1:atoms]
        eigenvalues = vector_template
        eigenvectors = matrix_template + I
        adiabatic_derivative = [matrix_template for _ in CartesianIndices(derivative)]
        nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(derivative)]
        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)

        new{T,M,n,n^2}(model,
            potential, derivative, eigenvalues, eigenvectors,
            adiabatic_derivative, nonadiabatic_coupling,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

mutable struct RingPolymerDiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Vector{Hermitian{T,SMatrix{S,S,T,L}}}
    derivative::Array{Hermitian{T,SMatrix{S,S,T,L}},3}
    eigenvalues::Vector{SVector{S,T}}
    eigenvectors::Vector{SMatrix{S,S,T,L}}
    adiabatic_derivative::Array{SMatrix{S,S,T,L},3}
    nonadiabatic_coupling::Array{SMatrix{S,S,T,L},3}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function RingPolymerDiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer, beads::Integer) where {T,M<:Model}
        n = Int(model.n_states)
        matrix_template = NonadiabaticModels.matrix_template(model, T)
        vector_template = NonadiabaticModels.vector_template(model, T)

        potential = [Hermitian(matrix_template) for _=1:beads]
        derivative = [Hermitian(matrix_template) for _=1:DoFs, _=1:atoms, _=1:beads]
        eigenvalues = [vector_template for _=1:beads]
        eigenvectors = [matrix_template + I for _=1:beads]
        adiabatic_derivative = [matrix_template for _=1:DoFs, _=1:atoms, _=1:beads]
        nonadiabatic_coupling = [matrix_template for _=1:DoFs, _=1:atoms, _=1:beads]
        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)
        new{T,M,n,n^2}(model, potential, derivative, eigenvalues, eigenvectors, adiabatic_derivative, nonadiabatic_coupling,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    DiabaticCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    AdiabaticCalculator{T}(model, DoFs, atoms)
end
function Calculator(model::DiabaticModel, DoFs::Integer, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerDiabaticCalculator{T}(model, DoFs, atoms, beads)
end
function Calculator(model::AdiabaticModel, DoFs::Integer, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerAdiabaticCalculator{T}(model, DoFs, atoms, beads)
end

evaluate_potential!(calc::AbstractCalculator, R) = calc.potential = potential(calc.model, R)

function evaluate_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        calc.potential[i] = potential(calc.model, R[:,:,i])
    end
end

function evaluate_centroid_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    calc.potential[1] = potential(calc.model, get_centroid(R))
end

evaluate_derivative!(calc::AbstractCalculator, R) = derivative!(calc.model, calc.derivative, R)

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        derivative!(calc.model, calc.derivative[:,:,i], R[:,:,i])
    end
end

function evaluate_centroid_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views derivative!(calc.model, calc.derivative[:,:,1], get_centroid(R))
end

function eigen!(calc::DiabaticCalculator)
    eig = eigen(calc.potential)
    calc.eigenvalues = eig.values
    calc.eigenvectors = correct_phase(eig.vectors, calc.eigenvectors)
    return nothing
end

function eigen!(calc::RingPolymerDiabaticCalculator)
    for i=1:length(calc.potential)
        eig = eigen(calc.potential[i])
        calc.eigenvalues[i] = eig.values
        calc.eigenvectors[i] = correct_phase(eig.vectors, calc.eigenvectors[i])
    end
    return nothing
end

function correct_phase(new_vectors::SMatrix, old_vectors::SMatrix)
    n = size(new_vectors, 1)
    vect = SVector{n}(sign(dot(new_vectors[:,i], old_vectors[:,i])) for i=1:n)
    return new_vectors .* vect'
end

function transform_derivative!(calc::AbstractDiabaticCalculator)
    for I in eachindex(calc.derivative)
        calc.adiabatic_derivative[I] = calc.eigenvectors' * calc.derivative[I] * calc.eigenvectors
    end
end

function transform_derivative!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.derivative, 3) # Beads
        for j in axes(calc.derivative, 2) # Atoms
            for k in axes(calc.derivative, 1) # DoFs
                calc.adiabatic_derivative[k,j,i] = calc.eigenvectors[i]' * calc.derivative[k,j,i] * calc.eigenvectors[i]
            end
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::AbstractDiabaticCalculator)
    for I in eachindex(calc.adiabatic_derivative)
        calc.nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(calc.adiabatic_derivative[I], calc.eigenvalues)
    end
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator)
    @views for i in axes(calc.nonadiabatic_coupling, 3) # Beads
        for I in CartesianIndices(size(calc.adiabatic_derivative)[1:2])
            calc.nonadiabatic_coupling[I,i] = evaluate_nonadiabatic_coupling(calc.adiabatic_derivative[I,i], calc.eigenvalues[i])
        end
    end
end

"""
# References

- HammesSchifferTully_JChemPhys_101_4657_1994 Eq. (32)
- SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, section 2.3
"""
function evaluate_nonadiabatic_coupling(adiabatic_derivative::SMatrix, eigenvalues::SVector)
    n = length(eigenvalues)
    SMatrix{n,n}(
        (i != j ? adiabatic_derivative[j,i] / (eigenvalues[i] - eigenvalues[j]) : 0
        for j=1:n, i=1:n))
end

function update_electronics!(calculator::AbstractDiabaticCalculator, r::AbstractArray)
    # nuclear DoF
    evaluate_potential!(calculator, r)
    # nuclear DoF
    evaluate_derivative!(calculator, r)
    # nuclear DoF; gets eigenvector
    eigen!(calculator)
    # nuclear DoF; calculates adiabatic derivative, bcs H set up diabatically:
    # eigenvect'*derivatives*eigenvect
    transform_derivative!(calculator)
    # electronic DoF (but based on nuclear forces and eigenvalues)
    evaluate_nonadiabatic_coupling!(calculator)
end

include("large_diabatic.jl")
include("friction.jl")

end # module