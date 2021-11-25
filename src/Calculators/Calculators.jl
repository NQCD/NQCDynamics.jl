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

using LinearAlgebra: LinearAlgebra, Hermitian, I, Eigen, tr
using StaticArrays: SMatrix, SVector

using NonadiabaticModels: NonadiabaticModels, Model, nstates
using NonadiabaticModels.AdiabaticModels: AdiabaticModel
using NonadiabaticModels.DiabaticModels: DiabaticModel, DiabaticFrictionModel
using NonadiabaticModels.FrictionModels: AdiabaticFrictionModel

using NonadiabaticMolecularDynamics: RingPolymers, ndofs

"""
    AbstractCalculator{M<:Model}

Top-level type for all calculators.

Each concrete calculator contains the `Model` and the fields to store the quantities
obtained from the model.
"""
abstract type AbstractCalculator{T,M<:Model} end
abstract type AbstractAdiabaticCalculator{T,M<:AdiabaticModel} <: AbstractCalculator{T,M} end
abstract type AbstractDiabaticCalculator{T,M<:Union{DiabaticFrictionModel,DiabaticModel}} <: AbstractCalculator{T,M} end
abstract type AbstractStaticDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{T,M} end
abstract type AbstractFrictionCalculator{T,M<:AdiabaticFrictionModel} <: AbstractCalculator{T,M} end

Base.eltype(::AbstractCalculator{T}) where {T} = T

mutable struct AdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{T,M}
    model::M
    potential::T
    derivative::Matrix{T}
    function AdiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, 0, zeros(ndofs(model), atoms))
    end
end

struct RingPolymerAdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{T,M}
    model::M
    potential::Vector{T}
    derivative::Array{T,3}
    function RingPolymerAdiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        new{T,M}(model, zeros(beads), zeros(ndofs(model), atoms, beads))
    end
end

mutable struct DiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::Hermitian{T,SMatrix{S,S,T,L}}
    derivative::Matrix{Hermitian{T,SMatrix{S,S,T,L}}}
    eigen::LinearAlgebra.Eigen{T,T,SMatrix{S,S,T,L},SVector{S,T}}
    adiabatic_derivative::Matrix{SMatrix{S,S,T,L}}
    nonadiabatic_coupling::Matrix{SMatrix{S,S,T,L}}
    tmp_position::Matrix{T}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function DiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        matrix_template = NonadiabaticModels.DiabaticModels.matrix_template(model, T)
        vector_template = NonadiabaticModels.DiabaticModels.vector_template(model, T)

        potential = Hermitian(matrix_template)
        derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms]
        eigen = Eigen(vector_template, matrix_template + I)
        adiabatic_derivative = [matrix_template for _ in CartesianIndices(derivative)]
        nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(derivative)]
        tmp_position = zeros(T, ndofs(model), atoms)
        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)

        new{T,M,n,n^2}(model,
            potential, derivative, eigen,
            adiabatic_derivative, nonadiabatic_coupling,
            tmp_position, tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

mutable struct RingPolymerDiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::Vector{Hermitian{T,SMatrix{S,S,T,L}}}
    traceless_potential::Vector{Hermitian{T,SMatrix{S,S,T,L}}}
    V̄::Vector{T}
    derivative::Array{Hermitian{T,SMatrix{S,S,T,L}},3}
    traceless_derivative::Array{Hermitian{T,SMatrix{S,S,T,L}},3}
    D̄::Array{T,3}
    eigen::Vector{LinearAlgebra.Eigen{T,T,SMatrix{S,S,T,L},SVector{S,T}}}
    adiabatic_derivative::Array{SMatrix{S,S,T,L},3}
    traceless_adiabatic_derivative::Array{SMatrix{S,S,T,L},3}
    nonadiabatic_coupling::Array{SMatrix{S,S,T,L},3}

    centroid_potential::Hermitian{T,SMatrix{S,S,T,L}}
    centroid_derivative::Matrix{Hermitian{T,SMatrix{S,S,T,L}}}
    centroid_eigen::LinearAlgebra.Eigen{T,T,SMatrix{S,S,T,L},SVector{S,T}}
    centroid_adiabatic_derivative::Matrix{SMatrix{S,S,T,L}}
    centroid_nonadiabatic_coupling::Matrix{SMatrix{S,S,T,L}}
    tmp_centroid::Matrix{T}

    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function RingPolymerDiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        n = nstates(model)
        matrix_template = NonadiabaticModels.DiabaticModels.matrix_template(model, T)
        vector_template = NonadiabaticModels.DiabaticModels.vector_template(model, T)

        potential = [Hermitian(matrix_template) for _=1:beads]
        traceless_potential = [Hermitian(matrix_template) for _=1:beads]
        V̄ = zeros(beads)

        derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        traceless_derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        D̄ = zeros(ndofs(model), atoms, beads)

        adiabatic_derivative = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]
        traceless_adiabatic_derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms, _=1:beads]

        eigen = [Eigen(vector_template, matrix_template + I) for _=1:beads]
        nonadiabatic_coupling = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]

        centroid_potential = Hermitian(matrix_template)
        centroid_derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms]
        centroid_eigen = Eigen(vector_template, matrix_template + I)
        centroid_adiabatic_derivative = [matrix_template for _ in CartesianIndices(centroid_derivative)]
        centroid_nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(centroid_derivative)]
        tmp_centroid = zeros(T, ndofs(model), atoms)

        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)
        new{T,M,n,n^2}(model,
            potential, traceless_potential, V̄,
            derivative, traceless_derivative, D̄,
            eigen, adiabatic_derivative, traceless_adiabatic_derivative, nonadiabatic_coupling,
            centroid_potential, centroid_derivative, centroid_eigen,
            centroid_adiabatic_derivative, centroid_nonadiabatic_coupling, tmp_centroid,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

function Calculator(model::DiabaticModel, atoms::Integer, t::Type{T}) where {T}
    DiabaticCalculator{t}(model, atoms)
end
function Calculator(model::AdiabaticModel, atoms::Integer, t::Type{T}) where {T}
    AdiabaticCalculator{t}(model, atoms)
end
function Calculator(model::DiabaticModel, atoms::Integer, beads::Integer, t::Type{T}) where {T}
    RingPolymerDiabaticCalculator{t}(model, atoms, beads)
end
function Calculator(model::AdiabaticModel, atoms::Integer, beads::Integer, t::Type{T}) where {T}
    RingPolymerAdiabaticCalculator{t}(model, atoms, beads)
end

function evaluate_potential!(calc::AbstractCalculator, R)
    calc.potential = NonadiabaticModels.potential(calc.model, R)
end

function evaluate_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        calc.potential[i] = NonadiabaticModels.potential(calc.model, R[:,:,i])
    end
end

function evaluate_V̄!(calc::RingPolymerDiabaticCalculator)
    for i in 1:length(calc.V̄)
        calc.V̄[i] = tr(calc.potential[i]) / nstates(calc.model)
    end
end

function evaluate_traceless_potential!(calc::RingPolymerDiabaticCalculator)
    n = nstates(calc.model)
    for I in 1:length(calc.traceless_potential)
        calc.traceless_potential[I] = Hermitian(SMatrix{n,n}(
            i != j ? calc.potential[I][j,i] : calc.potential[I][j,i] - calc.V̄[I] for j=1:n, i=1:n
        ))
    end
end

function evaluate_centroid_potential!(calc::AbstractCalculator, R::AbstractMatrix)
    calc.centroid_potential = NonadiabaticModels.potential(calc.model, R)
end

function evaluate_derivative!(calc::AbstractCalculator, R)
    NonadiabaticModels.derivative!(calc.model, calc.derivative, R)
end

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NonadiabaticModels.derivative!(calc.model, calc.derivative[:,:,i], R[:,:,i])
    end
end

function evaluate_D̄!(calc::RingPolymerDiabaticCalculator)
    for I in eachindex(calc.D̄)
        calc.D̄[I] = tr(calc.derivative[I]) / nstates(calc.model)
    end
end

function evaluate_traceless_derivative!(calc::RingPolymerDiabaticCalculator)
    n = nstates(calc.model)
    for I in eachindex(calc.derivative)
        calc.traceless_derivative[I] = Hermitian(SMatrix{n,n}(
            i != j ? calc.derivative[I][j,i] : calc.derivative[I][j,i] - calc.D̄[I] for j=1:n, i=1:n
        ))
    end
end

function evaluate_traceless_adiabatic_derivative!(calc::RingPolymerDiabaticCalculator)
    n = nstates(calc.model)
    for I in eachindex(calc.adiabatic_derivative)
        calc.traceless_adiabatic_derivative[I] = SMatrix{n,n}(
            i != j ? calc.adiabatic_derivative[I][j,i] : calc.adiabatic_derivative[I][j,i] - calc.D̄[I] for j=1:n, i=1:n
        )
    end
end

function evaluate_centroid_derivative!(calc::AbstractCalculator, R::AbstractMatrix)
    NonadiabaticModels.derivative!(calc.model, calc.centroid_derivative, R)
end

function eigen!(calc::DiabaticCalculator)
    eig = LinearAlgebra.eigen(calc.potential)
    corrected_vectors = correct_phase(eig.vectors, calc.eigen.vectors)
    calc.eigen = Eigen(eig.values, corrected_vectors)
    return nothing
end

function centroid_eigen!(calc::RingPolymerDiabaticCalculator)
    eig = LinearAlgebra.eigen(calc.centroid_potential)
    corrected_vectors = correct_phase(eig.vectors, calc.centroid_eigen.vectors)
    calc.centroid_eigen = Eigen(eig.values, corrected_vectors)
    return nothing
end

function eigen!(calc::RingPolymerDiabaticCalculator)
    for i=1:length(calc.potential)
        eig = LinearAlgebra.eigen(calc.potential[i])
        corrected_vectors = correct_phase(eig.vectors, calc.eigen[i].vectors)
        calc.eigen[i] = Eigen(eig.values, corrected_vectors)
    end
    return nothing
end

function correct_phase(new_vectors::SMatrix, old_vectors::SMatrix)
    n = size(new_vectors, 1)
    vect = SVector{n}(sign(LinearAlgebra.dot(new_vectors[:,i], old_vectors[:,i])) for i=1:n)
    return new_vectors .* vect'
end

function transform_derivative!(calc::AbstractDiabaticCalculator)
    for I in eachindex(calc.derivative)
        calc.adiabatic_derivative[I] = calc.eigen.vectors' * calc.derivative[I] * calc.eigen.vectors
    end
end

function transform_centroid_derivative!(calc::RingPolymerDiabaticCalculator)
    for I in eachindex(calc.centroid_derivative)
        calc.centroid_adiabatic_derivative[I] = calc.centroid_eigen.vectors' * calc.centroid_derivative[I] * calc.centroid_eigen.vectors
    end
end

function transform_derivative!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.derivative, 3) # Beads
        for j in axes(calc.derivative, 2) # Atoms
            for k in axes(calc.derivative, 1) # DoFs
                calc.adiabatic_derivative[k,j,i] = calc.eigen[i].vectors' * calc.derivative[k,j,i] * calc.eigen[i].vectors
            end
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::AbstractDiabaticCalculator)
    for I in eachindex(calc.adiabatic_derivative)
        calc.nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(calc.adiabatic_derivative[I], calc.eigen.values)
    end
end

function evaluate_centroid_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator)
    for I in eachindex(calc.centroid_adiabatic_derivative)
        calc.centroid_nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(calc.centroid_adiabatic_derivative[I], calc.centroid_eigen.values)
    end
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator)
    for i in axes(calc.nonadiabatic_coupling, 3) # Beads
        for I in CartesianIndices(size(calc.adiabatic_derivative)[1:2])
            calc.nonadiabatic_coupling[I,i] = evaluate_nonadiabatic_coupling(calc.adiabatic_derivative[I,i], calc.eigen[i].values)
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

"""
Evaluates all electronic properties for the current position `r`.

# Properties evaluated:
- Diabatic potential
- Diabatic derivative
- Eigenvalues and eigenvectors
- Adiabatic derivative
- Nonadiabatic coupling
"""
function update_electronics!(calculator::AbstractDiabaticCalculator, r::AbstractArray)
    evaluate_potential!(calculator, r)
    evaluate_derivative!(calculator, r)
    eigen!(calculator)
    transform_derivative!(calculator)
    evaluate_nonadiabatic_coupling!(calculator)
end

function update_electronics!(calculator::RingPolymerDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    evaluate_potential!(calculator, r)
    evaluate_derivative!(calculator, r)
    eigen!(calculator)
    transform_derivative!(calculator)
    evaluate_nonadiabatic_coupling!(calculator)

    update_centroid_electronics!(calculator, r)
end

function update_centroid_electronics!(calculator::RingPolymerDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    RingPolymers.get_centroid!(calculator.tmp_centroid, r)
    evaluate_centroid_potential!(calculator, calculator.tmp_centroid)
    evaluate_centroid_derivative!(calculator, calculator.tmp_centroid)
    centroid_eigen!(calculator)
    transform_centroid_derivative!(calculator)
    evaluate_centroid_nonadiabatic_coupling!(calculator)
end

include("large_diabatic.jl")
include("friction.jl")

position(calc) = calc.tmp_position

for field ∈ (:potential, :derivative)
    get_field = Symbol(:get, :_, field)
    evaluate_field = Symbol(:evaluate, :_, field, :!)

    eval(quote
        function $get_field(calc, r)
            if position(calc) ≈ r
                return calc.$field
            else
                copyto!(position(calc), r)
                return $evaluate_field(calc, r)
            end
        end
    end)
end

end # module
