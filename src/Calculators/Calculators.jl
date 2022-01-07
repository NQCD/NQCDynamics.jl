"""
    Calculators

This module exists to bridge the gap between the `Models` and the `Dynamics`.

Here, we provide functions and types for evaluating and storing quantities obtained from the
`Models`.
In addition any further manipulation of those quantities, such as computing eigenvalues,
is included here.

This module is largely needed to facilitate integration of both ring polymer and
classical dynamics to allow using the same models and functions for both.
Specific ring polymer types are provided that have the extra fields and methods needed
to evaluate the quantities for each bead. 
"""
module Calculators

using LinearAlgebra: LinearAlgebra, Hermitian, I, Eigen, tr
using StaticArrays: SMatrix, SVector

using NQCModels: NQCModels, Model, nstates
using NQCModels.AdiabaticModels: AdiabaticModel
using NQCModels.DiabaticModels: DiabaticModel, DiabaticFrictionModel
using NQCModels.FrictionModels: AdiabaticFrictionModel

using NQCDynamics: RingPolymers, ndofs

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

mutable struct DependentField{V,R}
    value::V
    position::R
end
needsupdate(field::DependentField, r) = field.position != r

function Base.getproperty(calc::AbstractCalculator, name::Symbol)
    field = getfield(calc, name)
    if field isa DependentField
        return field.value
    else
        return field
    end
end

function Base.setproperty!(calc::AbstractCalculator, name::Symbol, x)
    field = getfield(calc, name)
    if field isa DependentField
        setfield!(field, :value, x)
    else
        setfield!(calc, name, x)
    end
end

"""
Each of the quantities specified here has functions:
`get_quantity(calculator, r)`
`evaluate_quantity(calculator, r)!`

The user should access only the former.
This will ensure quantities are correctly evaluated and cached accordingly.

The latter is called by the former and is where the details required to calculate the quantity are found.
"""
const quantities = [
    :potential,
    :derivative,
    :eigen,
    :adiabatic_derivative,
    :nonadiabatic_coupling,

    :traceless_potential,
    :V̄,
    :traceless_derivative,
    :D̄,
    :traceless_adiabatic_derivative,

    :centroid,
    :centroid_potential,
    :centroid_derivative,
    :centroid_eigen,
    :centroid_adiabatic_derivative,
    :centroid_nonadiabatic_coupling,

    :friction,
]

for quantity in quantities
    get_quantity = Symbol(:get_, quantity)
    evaluate_quantity! = Symbol(:evaluate_, quantity, :!)
    field = Expr(:call, :getfield, :calculator, QuoteNode(quantity))

    @eval function $(get_quantity)(calculator, r)

        if needsupdate($field, r)
            copyto!($field.position, r)
            $(evaluate_quantity!)(calculator, r)
        end
        return $field.value
    end
end

struct AdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{T,M}
    model::M
    potential::DependentField{T,Matrix{T}}
    derivative::DependentField{Matrix{T},Matrix{T}}
    function AdiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        potential = zero(T)
        derivative = zeros(T, ndofs(model), atoms)
        position = fill(NaN, ndofs(model), atoms)
        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position))
        )
    end
end

struct RingPolymerAdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{T,M}
    model::M
    potential::DependentField{Vector{T},Array{T,3}}
    derivative::DependentField{Array{T,3},Array{T,3}}
    function RingPolymerAdiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        potential = zeros(T, beads)
        derivative = zeros(T, ndofs(model), atoms, beads)
        position = fill(NaN, ndofs(model), atoms, beads)
        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position))
        )
    end
end

struct DiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::DependentField{Hermitian{T,SMatrix{S,S,T,L}},Matrix{T}}
    derivative::DependentField{Matrix{Hermitian{T,SMatrix{S,S,T,L}}},Matrix{T}}
    eigen::DependentField{LinearAlgebra.Eigen{T,T,SMatrix{S,S,T,L},SVector{S,T}},Matrix{T}}
    adiabatic_derivative::DependentField{Matrix{SMatrix{S,S,T,L}},Matrix{T}}
    nonadiabatic_coupling::DependentField{Matrix{SMatrix{S,S,T,L}},Matrix{T}}
    function DiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        matrix_template = NQCModels.DiabaticModels.matrix_template(model, T)
        vector_template = NQCModels.DiabaticModels.vector_template(model, T)

        potential = Hermitian(matrix_template)
        derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms]
        eigen = Eigen(vector_template, matrix_template + I)
        adiabatic_derivative = [matrix_template for _ in CartesianIndices(derivative)]
        nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(derivative)]
        position = fill(NaN, ndofs(model), atoms)

        new{T,M,n,n^2}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),
        )
    end
end

struct RingPolymerDiabaticCalculator{T,M,S,L} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::DependentField{Vector{Hermitian{T,SMatrix{S,S,T,L}}},Array{T,3}}
    derivative::DependentField{Array{Hermitian{T,SMatrix{S,S,T,L}},3},Array{T,3}}
    eigen::DependentField{Vector{LinearAlgebra.Eigen{T,T,SMatrix{S,S,T,L},SVector{S,T}}},Array{T,3}}
    adiabatic_derivative::DependentField{Array{SMatrix{S,S,T,L},3},Array{T,3}}
    nonadiabatic_coupling::DependentField{Array{SMatrix{S,S,T,L},3},Array{T,3}}

    traceless_potential::DependentField{Vector{Hermitian{T,SMatrix{S,S,T,L}}},Array{T,3}}
    V̄::DependentField{Vector{T},Array{T,3}}
    traceless_derivative::DependentField{Array{Hermitian{T,SMatrix{S,S,T,L}},3},Array{T,3}}
    D̄::DependentField{Array{T,3},Array{T,3}}
    traceless_adiabatic_derivative::DependentField{Array{SMatrix{S,S,T,L},3},Array{T,3}}

    centroid::DependentField{Matrix{T},Array{T,3}}
    centroid_potential::DependentField{Hermitian{T,SMatrix{S,S,T,L}},Array{T,3}}
    centroid_derivative::DependentField{Matrix{Hermitian{T,SMatrix{S,S,T,L}}},Array{T,3}}
    centroid_eigen::DependentField{LinearAlgebra.Eigen{T,T,SMatrix{S,S,T,L},SVector{S,T}},Array{T,3}}
    centroid_adiabatic_derivative::DependentField{Matrix{SMatrix{S,S,T,L}},Array{T,3}}
    centroid_nonadiabatic_coupling::DependentField{Matrix{SMatrix{S,S,T,L}},Array{T,3}}
    function RingPolymerDiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        n = nstates(model)
        matrix_template = NQCModels.DiabaticModels.matrix_template(model, T)
        vector_template = NQCModels.DiabaticModels.vector_template(model, T)

        potential = [Hermitian(matrix_template) for _=1:beads]
        derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        adiabatic_derivative = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]
        eigen = [Eigen(vector_template, matrix_template + I) for _=1:beads]
        nonadiabatic_coupling = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]

        traceless_potential = [Hermitian(matrix_template) for _=1:beads]
        V̄ = zeros(beads)
        traceless_derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        D̄ = zeros(ndofs(model), atoms, beads)
        traceless_adiabatic_derivative = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]

        centroid = zeros(T, ndofs(model), atoms)
        centroid_potential = Hermitian(matrix_template)
        centroid_derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms]
        centroid_eigen = Eigen(vector_template, matrix_template + I)
        centroid_adiabatic_derivative = [matrix_template for _ in CartesianIndices(centroid_derivative)]
        centroid_nonadiabatic_coupling = [matrix_template for _ in CartesianIndices(centroid_derivative)]

        position = fill(NaN, ndofs(model), atoms, beads)

        new{T,M,n,n^2}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),

            DependentField(traceless_potential, copy(position)),
            DependentField(V̄, copy(position)),
            DependentField(traceless_derivative, copy(position)),
            DependentField(D̄, copy(position)),
            DependentField(traceless_adiabatic_derivative, copy(position)),

            DependentField(centroid, copy(position)),
            DependentField(centroid_potential, copy(position)),
            DependentField(centroid_derivative, copy(position)),
            DependentField(centroid_eigen, copy(position)),
            DependentField(centroid_adiabatic_derivative, copy(position)),
            DependentField(centroid_nonadiabatic_coupling, copy(position)),
        )
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
    calc.potential = NQCModels.potential(calc.model, R)
end

function evaluate_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        calc.potential[i] = NQCModels.potential(calc.model, R[:,:,i])
    end
end

function evaluate_V̄!(calc::RingPolymerDiabaticCalculator, r)
    potential = get_potential(calc, r)
    for i in 1:length(calc.V̄)
        calc.V̄[i] = tr(potential[i]) / nstates(calc.model)
    end
end

function evaluate_traceless_potential!(calc::RingPolymerDiabaticCalculator, r)
    n = nstates(calc.model)
    potential = get_potential(calc, r)
    V̄ = get_V̄(calc, r)
    for I in eachindex(potential)
        calc.traceless_potential[I] = Hermitian(SMatrix{n,n}(
            i != j ? potential[I][j,i] : potential[I][j,i] - V̄[I] for j=1:n, i=1:n
        ))
    end
end

function evaluate_centroid!(calc::AbstractCalculator, r::Array{T,3}) where {T}
    RingPolymers.get_centroid!(calc.centroid, r)
end

function evaluate_centroid_potential!(calc::AbstractCalculator, r::Array{T,3}) where {T}
    centroid = get_centroid(calc, r)
    calc.centroid_potential = NQCModels.potential(calc.model, centroid)
end

function evaluate_derivative!(calc::AbstractCalculator, R)
    NQCModels.derivative!(calc.model, calc.derivative, R)
end

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NQCModels.derivative!(calc.model, calc.derivative[:,:,i], R[:,:,i])
    end
end

function evaluate_D̄!(calc::RingPolymerDiabaticCalculator, r)
    derivative = get_derivative(calc, r)
    for I in eachindex(derivative)
        calc.D̄[I] = tr(derivative[I]) / nstates(calc.model)
    end
end

function evaluate_traceless_derivative!(calc::RingPolymerDiabaticCalculator, r)
    n = nstates(calc.model)
    derivative = get_derivative(calc, r)
    D̄ = get_D̄(calc, r)
    for I in eachindex(derivative)
        calc.traceless_derivative[I] = Hermitian(SMatrix{n,n}(
            i != j ? derivative[I][j,i] : derivative[I][j,i] - D̄[I] for j=1:n, i=1:n
        ))
    end
end

function evaluate_traceless_adiabatic_derivative!(calc::RingPolymerDiabaticCalculator, r)
    n = nstates(calc.model)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)
    D̄ = get_D̄(calc, r)
    for I in eachindex(D̄)
        calc.traceless_adiabatic_derivative[I] = SMatrix{n,n}(
            i != j ? adiabatic_derivative[I][j,i] : adiabatic_derivative[I][j,i] - D̄[I] for j=1:n, i=1:n
        )
    end
end

function evaluate_centroid_derivative!(calc::AbstractCalculator, r::AbstractArray{T,3}) where {T}
    centroid = get_centroid(calc, r)
    NQCModels.derivative!(calc.model, calc.centroid_derivative, centroid)
end

function evaluate_eigen!(calc::AbstractDiabaticCalculator, r)
    potential = get_potential(calc, r)
    eig = LinearAlgebra.eigen(potential)
    corrected_vectors = correct_phase(eig.vectors, calc.eigen.vectors)
    calc.eigen = Eigen(eig.values, corrected_vectors)
    return nothing
end

function evaluate_centroid_eigen!(calc::RingPolymerDiabaticCalculator, r)
    potential = get_centroid_potential(calc, r)
    eig = LinearAlgebra.eigen(potential)
    corrected_vectors = correct_phase(eig.vectors, calc.centroid_eigen.vectors)
    calc.centroid_eigen = Eigen(eig.values, corrected_vectors)
    return nothing
end

function evaluate_eigen!(calc::RingPolymerDiabaticCalculator, r)
    potential = get_potential(calc, r)
    for i=1:length(potential)
        eig = LinearAlgebra.eigen(potential[i])
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

function evaluate_adiabatic_derivative!(calc::AbstractDiabaticCalculator, r)
    U = get_eigen(calc, r).vectors
    diabatic_derivative = get_derivative(calc, r)
    for I in eachindex(diabatic_derivative)
        calc.adiabatic_derivative[I] = U' * diabatic_derivative[I] * U
    end
end

function evaluate_centroid_adiabatic_derivative!(calc::RingPolymerDiabaticCalculator, r)
    centroid_derivative = get_centroid_derivative(calc, r)
    centroid_eigen = get_centroid_eigen(calc, r)
    for I in eachindex(centroid_derivative)
        calc.centroid_adiabatic_derivative[I] = centroid_eigen.vectors' * centroid_derivative[I] * centroid_eigen.vectors
    end
end

function evaluate_adiabatic_derivative!(calc::RingPolymerDiabaticCalculator, r)
    derivative = get_derivative(calc, r)
    eigen = get_eigen(calc, r)
    for i in axes(derivative, 3) # Beads
        for j in axes(derivative, 2) # Atoms
            for k in axes(derivative, 1) # DoFs
                calc.adiabatic_derivative[k,j,i] = eigen[i].vectors' * derivative[k,j,i] * eigen[i].vectors
            end
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::AbstractDiabaticCalculator, r)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)
    eigen = get_eigen(calc, r)
    for I in eachindex(calc.adiabatic_derivative)
        calc.nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(adiabatic_derivative[I], eigen.values)
    end
end

function evaluate_centroid_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator, r)
    centroid_adiabatic_derivative = get_centroid_adiabatic_derivative(calc, r)
    centroid_eigen = get_centroid_eigen(calc, r)
    for I in eachindex(centroid_adiabatic_derivative)
        calc.centroid_nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(centroid_adiabatic_derivative[I], centroid_eigen.values)
    end
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator, r)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)
    eigen = get_eigen(calc, r)
    for i in eachindex(eigen) # Beads
        for I in CartesianIndices(size(adiabatic_derivative)[1:2])
            calc.nonadiabatic_coupling[I,i] = evaluate_nonadiabatic_coupling(adiabatic_derivative[I,i], eigen[i].values)
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

include("friction.jl")
include("large_diabatic.jl")

end # module
