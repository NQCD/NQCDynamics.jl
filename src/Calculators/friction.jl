
mutable struct FrictionCalculator{T,M} <: AbstractFrictionCalculator{T,M}
    model::M
    potential::DependentField{T,Matrix{T}}
    derivative::DependentField{Matrix{T},Matrix{T}}
    friction::DependentField{Matrix{T},Matrix{T}}
    stats::Dict{Symbol,Int}
    function FrictionCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}

        potential = zero(T)
        derivative = zeros(ndofs(model), atoms)
        friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms)

        position = fill(NaN, ndofs(model), atoms)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
            :friction=>0,
        )

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(friction, copy(position)),
            stats
        )
    end
end

struct RingPolymerFrictionCalculator{T,M} <: AbstractFrictionCalculator{T,M}
    model::M
    potential::DependentField{Vector{T},Array{T,3}}
    derivative::DependentField{Array{T,3},Array{T,3}}
    friction::DependentField{Array{T,3},Array{T,3}}
    stats::Dict{Symbol,Int}
    function RingPolymerFrictionCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}

        potential = zeros(beads)
        derivative = zeros(ndofs(model), atoms, beads)
        friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms, beads)

        position = fill(NaN, ndofs(model), atoms, beads)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
            :friction=>0,
        )

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(friction, copy(position)),
            stats
        )
    end
end

function Calculator(model::AdiabaticFrictionModel, atoms::Integer, t::Type{T}) where {T}
    FrictionCalculator{t}(model, atoms)
end
function Calculator(model::AdiabaticFrictionModel, atoms::Integer, beads::Integer, t::Type{T}) where {T}
    RingPolymerFrictionCalculator{t}(model, atoms, beads)
end
function Calculator(model::CompositeModel, atoms::Integer, t::Type{T}) where {T}
    if any([!isa(s.model, AdiabaticModel) for s in NQCModels.get_pes_models(model.subsystems)])
        throw(ArgumentError("Currently, only CompositeModels using AdiabaticModels to supply a PES are supported. "))
    end
    FrictionCalculator{t}(model, atoms)
end

function Calculator(model::CompositeModel, atoms::Integer, beads::Integer, t::Type{T}) where {T}
    if any([!isa(s.model, AdiabaticModel) for s in NQCModels.get_pes_models(model.subsystems)])
        throw(ArgumentError("Currently, only CompositeModels using AdiabaticModels to supply a PES are supported. "))
    end
    RingPolymerFrictionCalculator{t}(model, atoms, beads)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractMatrix)
    calc.stats[:friction] += 1
    NQCModels.friction!(calc.model, calc.friction, R)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractArray{T,3}) where {T}
    calc.stats[:friction] += 1
    @views for i in axes(R, 3)
        NQCModels.friction!(calc.model, calc.friction[:,:,i], R[:,:,i])
    end
end
