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
using RingPolymerArrays: get_centroid!

using NQCModels: NQCModels, Model, nstates, mobileatoms, dofs, Subsystem, CompositeModel
using NQCModels.AdiabaticModels: AdiabaticModel
using NQCModels.DiabaticModels: DiabaticModel, DiabaticFrictionModel
using NQCModels.FrictionModels: AdiabaticFrictionModel

using NQCDynamics: ndofs

"""
    AbstractCalculator{M<:Model}

Top-level type for all calculators.

Each concrete calculator contains the `Model` and the fields to store the quantities
obtained from the model.
"""
abstract type AbstractCalculator{T,M<:Model} end
abstract type AbstractAdiabaticCalculator{T,M<:Union{AdiabaticModel, CompositeModel}} <: AbstractCalculator{T,M} end
abstract type AbstractDiabaticCalculator{T,M<:Union{DiabaticFrictionModel,DiabaticModel}} <: AbstractCalculator{T,M} end
abstract type AbstractStaticDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{T,M} end
abstract type AbstractFrictionCalculator{T,M<:Union{AdiabaticFrictionModel, CompositeModel}} <: AbstractCalculator{T,M} end

NQCModels.nstates(calc::AbstractCalculator) = NQCModels.nstates(calc.model)
NQCModels.eachstate(calc::AbstractCalculator) = NQCModels.eachstate(calc.model)
NQCModels.nelectrons(calc::AbstractCalculator) = NQCModels.nelectrons(calc.model)
NQCModels.eachelectron(calc::AbstractCalculator) = NQCModels.eachelectron(calc.model)
NQCModels.mobileatoms(calc::AbstractCalculator) = NQCModels.mobileatoms(calc.model, size(calc.derivative, 2))
NQCModels.dofs(calc::AbstractCalculator) = NQCModels.dofs(calc.model)
beads(calc) = Base.OneTo(length(calc.potential))
NQCModels.fermilevel(calc::AbstractCalculator) = NQCModels.fermilevel(calc.model)

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
`evaluate_quantity!(calculator, r)!`

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
    stats::Dict{Symbol,Int}
    function AdiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        potential = zero(T)
        derivative = zeros(T, ndofs(model), atoms)
        position = fill(NaN, ndofs(model), atoms)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
        )

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            stats
        )
    end
end

struct RingPolymerAdiabaticCalculator{T,M} <: AbstractAdiabaticCalculator{T,M}
    model::M
    potential::DependentField{Vector{T},Array{T,3}}
    derivative::DependentField{Array{T,3},Array{T,3}}
    stats::Dict{Symbol,Int}
    function RingPolymerAdiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        potential = zeros(T, beads)
        derivative = zeros(T, ndofs(model), atoms, beads)
        position = fill(NaN, ndofs(model), atoms, beads)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
        )

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            stats
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
    stats::Dict{Symbol,Int}
    function DiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        mat = NQCModels.DiabaticModels.matrix_template(model, T)
        vec = NQCModels.DiabaticModels.vector_template(model, T)

        potential = Hermitian(zero(mat))
        derivative = [Hermitian(zero(mat)) for _=1:ndofs(model), _=1:atoms]
        eigen = Eigen(zero(vec), zero(mat) + I)
        adiabatic_derivative = [zero(mat) for _ in CartesianIndices(derivative)]
        nonadiabatic_coupling = [zero(mat) for _ in CartesianIndices(derivative)]
        position = fill(NaN, ndofs(model), atoms)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
            :eigen=>0,
            :adiabatic_derivative=>0,
            :nonadiabatic_coupling=>0
        )

        new{T,M,n,n^2}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),
            stats
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

    stats::Dict{Symbol,Int}
    function RingPolymerDiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        n = nstates(model)
        mat = NQCModels.DiabaticModels.matrix_template(model, T)
        vec = NQCModels.DiabaticModels.vector_template(model, T)

        potential = [Hermitian(zero(mat)) for _=1:beads]
        derivative = [Hermitian(zero(mat)) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        adiabatic_derivative = [zero(mat) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        eigen = [Eigen(zero(vec), zero(mat) + I) for _=1:beads]
        nonadiabatic_coupling = [zero(mat) for _=1:ndofs(model), _=1:atoms, _=1:beads]

        traceless_potential = [Hermitian(zero(mat)) for _=1:beads]
        V̄ = zeros(beads)
        traceless_derivative = [Hermitian(zero(mat)) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        D̄ = zeros(ndofs(model), atoms, beads)
        traceless_adiabatic_derivative = [zero(mat) for _=1:ndofs(model), _=1:atoms, _=1:beads]

        centroid = zeros(T, ndofs(model), atoms)
        centroid_potential = Hermitian(zero(mat))
        centroid_derivative = [Hermitian(zero(mat)) for _=1:ndofs(model), _=1:atoms]
        centroid_eigen = Eigen(zero(vec), zero(mat) + I)
        centroid_adiabatic_derivative = [zero(mat) for _ in CartesianIndices(centroid_derivative)]
        centroid_nonadiabatic_coupling = [zero(mat) for _ in CartesianIndices(centroid_derivative)]

        position = fill(NaN, ndofs(model), atoms, beads)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
            :eigen=>0,
            :adiabatic_derivative=>0,
            :nonadiabatic_coupling=>0,
            :traceless_potential=>0,
            :V̄=>0,
            :traceless_derivative=>0,
            :D̄=>0,
            :traceless_adiabatic_derivative=>0,
            :centroid=>0,
            :centroid_potential=>0,
            :centroid_derivative=>0,
            :centroid_eigen=>0,
            :centroid_adiabatic_derivative=>0,
            :centroid_nonadiabatic_coupling=>0,
        )

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
            stats
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
    calc.stats[:potential] += 1
    calc.potential = NQCModels.potential(calc.model, R)
    return nothing
end

function evaluate_potential!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    calc.stats[:potential] += 1
    @views for i in axes(R, 3)
        calc.potential[i] = NQCModels.potential(calc.model, R[:,:,i])
    end
    return nothing
end

function evaluate_V̄!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:V̄] += 1
    potential = get_potential(calc, r)
    for i in 1:length(calc.V̄)
        calc.V̄[i] = tr(potential[i]) / nstates(calc.model)
    end
end

function evaluate_traceless_potential!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:traceless_potential] += 1
    n = nstates(calc.model)
    potential = get_potential(calc, r)
    V̄ = get_V̄(calc, r)
    for I in eachindex(potential)
        calc.traceless_potential[I] = Hermitian(SMatrix{n,n}(
            i != j ? potential[I][j,i] : potential[I][j,i] - V̄[I] for j=1:n, i=1:n
        ))
    end
end

function evaluate_centroid!(calc::AbstractCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:centroid] += 1
    get_centroid!(calc.centroid, r)
end

function evaluate_centroid_potential!(calc::AbstractCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:centroid_potential] += 1
    centroid = get_centroid(calc, r)
    calc.centroid_potential = NQCModels.potential(calc.model, centroid)
end

function evaluate_derivative!(calc::AbstractCalculator, R)
    calc.stats[:derivative] += 1
    NQCModels.derivative!(calc.model, calc.derivative, R)
end

function evaluate_derivative!(calc::AbstractCalculator, R::AbstractArray{T,3}) where {T}
    calc.stats[:derivative] += 1
    @views for i in axes(R, 3)
        NQCModels.derivative!(calc.model, calc.derivative[:,:,i], R[:,:,i])
    end
end

function evaluate_D̄!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:D̄] += 1
    derivative = get_derivative(calc, r)
    for I in eachindex(derivative)
        calc.D̄[I] = tr(derivative[I]) / nstates(calc.model)
    end
end

function evaluate_traceless_derivative!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:traceless_derivative] += 1
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
    calc.stats[:traceless_adiabatic_derivative] += 1
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
    calc.stats[:centroid_derivative] += 1
    centroid = get_centroid(calc, r)
    NQCModels.derivative!(calc.model, calc.centroid_derivative, centroid)
end

function evaluate_eigen!(calc::AbstractDiabaticCalculator, r)
    calc.stats[:eigen] += 1
    potential = get_potential(calc, r)
    eig = LinearAlgebra.eigen(potential)
    corrected_vectors = correct_phase(eig.vectors, calc.eigen.vectors)
    calc.eigen = Eigen(eig.values, corrected_vectors)
    return nothing
end

function evaluate_centroid_eigen!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:centroid_eigen] += 1
    potential = get_centroid_potential(calc, r)
    eig = LinearAlgebra.eigen(potential)
    corrected_vectors = correct_phase(eig.vectors, calc.centroid_eigen.vectors)
    calc.centroid_eigen = Eigen(eig.values, corrected_vectors)
    return nothing
end

function evaluate_eigen!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:eigen] += 1
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
    calc.stats[:adiabatic_derivative] += 1
    U = get_eigen(calc, r).vectors
    diabatic_derivative = get_derivative(calc, r)
    for I in eachindex(diabatic_derivative)
        calc.adiabatic_derivative[I] = U' * diabatic_derivative[I] * U
    end
end

function evaluate_centroid_adiabatic_derivative!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:centroid_adiabatic_derivative] += 1
    centroid_derivative = get_centroid_derivative(calc, r)
    centroid_eigen = get_centroid_eigen(calc, r)
    for I in eachindex(centroid_derivative)
        calc.centroid_adiabatic_derivative[I] = centroid_eigen.vectors' * centroid_derivative[I] * centroid_eigen.vectors
    end
end

function evaluate_adiabatic_derivative!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:adiabatic_derivative] += 1
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
    calc.stats[:nonadiabatic_coupling] += 1
    adiabatic_derivative = get_adiabatic_derivative(calc, r)
    eigen = get_eigen(calc, r)
    for I in eachindex(calc.adiabatic_derivative)
        calc.nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(adiabatic_derivative[I], eigen.values)
    end
end

function evaluate_centroid_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:centroid_nonadiabatic_coupling] += 1
    centroid_adiabatic_derivative = get_centroid_adiabatic_derivative(calc, r)
    centroid_eigen = get_centroid_eigen(calc, r)
    for I in eachindex(centroid_adiabatic_derivative)
        calc.centroid_nonadiabatic_coupling[I] = evaluate_nonadiabatic_coupling(centroid_adiabatic_derivative[I], centroid_eigen.values)
    end
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerDiabaticCalculator, r)
    calc.stats[:nonadiabatic_coupling] += 1
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
include("ring_polymer_large_diabatic.jl")

"""
Evaluates all electronic properties for the current position `r`.
# Properties evaluated:
- Diabatic potential
- Diabatic derivative
- Eigenvalues and eigenvectors
- Adiabatic derivative
- Nonadiabatic coupling

This should no longer be used, instead access the quantities directly with `get_quantity(calc, r)`.
"""
function update_electronics!(calculator::AbstractDiabaticCalculator, r::AbstractArray)
    get_nonadiabatic_coupling(calculator, r)
    return nothing
end

function update_electronics!(calculator::RingPolymerDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    get_nonadiabatic_coupling(calculator, r)
    update_centroid_electronics!(calculator, r)
    return nothing
end

function update_centroid_electronics!(calculator::RingPolymerDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    get_centroid_nonadiabatic_coupling(calculator, r)
    return nothing
end

end # module
