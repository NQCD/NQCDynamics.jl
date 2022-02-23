
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

# struct DiabaticFrictionCalculator{T,M} <: AbstractDiabaticCalculator{T,M}
#     model::M
#     potential::DependentField{Hermitian{T,Matrix{T}},Matrix{T}}
#     derivative::DependentField{Matrix{Hermitian{T,Matrix{T}}},Matrix{T}}
#     eigen::DependentField{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}},Matrix{T}}
#     adiabatic_derivative::DependentField{Matrix{Matrix{T}},Matrix{T}}
#     nonadiabatic_coupling::DependentField{Matrix{Matrix{T}},Matrix{T}}
#     friction::DependentField{Matrix{T},Matrix{T}}
#     tmp_mat::Matrix{T}

#     function DiabaticFrictionCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
#         mat = NQCModels.DiabaticModels.matrix_template(model, T)
#         vec = NQCModels.DiabaticModels.vector_template(model, T)

#         potential = Hermitian(zero(mat))
#         derivative = [Hermitian(zero(mat)) for i=1:ndofs(model), j=1:atoms]
#         eigen = Eigen(zero(vec), zero(mat)+I)
#         adiabatic_derivative = [zero(mat) for i=1:ndofs(model), j=1:atoms]
#         nonadiabatic_coupling = [zero(mat) for i=1:ndofs(model), j=1:atoms]
#         friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms)
#         tmp_mat = zero(mat)

#         position = fill(NaN, ndofs(model), atoms)

#         new{T,M}(
#             model,
#             DependentField(potential, copy(position)),
#             DependentField(derivative, copy(position)),
#             DependentField(eigen, copy(position)),
#             DependentField(adiabatic_derivative, copy(position)),
#             DependentField(nonadiabatic_coupling, copy(position)),
#             DependentField(friction, copy(position)),
#             tmp_mat
#         )
#     end
# end

# struct RingPolymerDiabaticFrictionCalculator{T,M} <: AbstractDiabaticCalculator{T,M}
#     model::M
#     potential::DependentField{Vector{Hermitian{T,Matrix{T}}},Array{T,3}}
#     derivative::DependentField{Array{Hermitian{T,Matrix{T}},3},Array{T,3}}
#     eigen::DependentField{Vector{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}}},Array{T,3}}
#     adiabatic_derivative::DependentField{Array{Matrix{T},3},Array{T,3}}
#     nonadiabatic_coupling::DependentField{Array{Matrix{T},3},Array{T,3}}
#     friction::DependentField{Matrix{T},Array{T,3}}
#     tmp_mat::Matrix{T}

#     function RingPolymerDiabaticFrictionCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
#         n = nstates(model)
#         matrix_template = NQCModels.DiabaticModels.matrix_template(model, T)
#         vector_template = NQCModels.DiabaticModels.vector_template(model, T)

#         potential = [Hermitian(matrix_template) for _=1:beads]
#         derivative = [Hermitian(matrix_template) for _=1:ndofs(model), _=1:atoms, _=1:beads]
#         adiabatic_derivative = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]
#         eigen = [Eigen(vector_template, matrix_template + I) for _=1:beads]
#         nonadiabatic_coupling = [matrix_template for _=1:ndofs(model), _=1:atoms, _=1:beads]
#         friction = zeros(ndofs(model)*atoms*beads, ndofs(model)*atoms*beads)
#         tmp_mat = zeros(T, n, n)

#         position = fill(NaN, ndofs(model), atoms, beads)

#         new{T,M}(
#             model,
#             DependentField(potential, copy(position)),
#             DependentField(derivative, copy(position)),
#             DependentField(eigen, copy(position)),
#             DependentField(adiabatic_derivative, copy(position)),
#             DependentField(nonadiabatic_coupling, copy(position)),
#             DependentField(friction, copy(position)),
#             tmp_mat
#         )
#     end
# end

function Calculator(model::AdiabaticFrictionModel, atoms::Integer, t::Type{T}) where {T}
    FrictionCalculator{t}(model, atoms)
end
# function Calculator(model::DiabaticFrictionModel, atoms::Integer, t::Type{T}) where {T}
#     DiabaticFrictionCalculator{t}(model, atoms)
# end
# function Calculator(model::DiabaticFrictionModel, atoms::Integer, beads::Integer, t::Type{T}) where {T}
#     RingPolymerDiabaticFrictionCalculator{t}(model, atoms, beads)
# end
function Calculator(model::AdiabaticFrictionModel, atoms::Integer, beads::Integer, t::Type{T}) where {T}
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

# @doc raw"""
#     evaluate_friction!(calc::DiabaticFrictionCalculator, R::AbstractMatrix)

# Evaluate the electronic friction for a model given in the diabatic representation.

# ```math
# γ = 2πħ ∑ᵢⱼ <i|dH|j><j|dH|i> (f(ϵᵢ)-f(ϵⱼ)) δ(ϵⱼ-ϵᵢ) / (ϵⱼ-ϵᵢ)
# ```
# Note that the delta function is approximated by a normalised gaussian.
# """
# function evaluate_friction!(calc::DiabaticFrictionCalculator, r::AbstractMatrix, β, σ, deriv=false)

#     fermi(ϵ, μ, β) = 1 / (1 + exp(β*(ϵ-μ)))
#     ∂fermi(ϵ, μ, β) = -β * exp(β*(ϵ-μ)) / (1 + exp(β*(ϵ-μ)))^2
#     gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))

#     adiabatic_derivative = get_adiabatic_derivative(calc, r)
#     eigen = get_eigen(calc, r)
#     μ = NQCModels.fermilevel(calc.model)

#     fill!(calc.friction, zero(eltype(calc)))
#     for I in eachindex(r)
#         for J in eachindex(r)
#             for m=1:nstates(calc.model)
#                 for n=m+1:nstates(calc.model)
#                     ϵₙ = eigen.values[n]
#                     ϵₘ = eigen.values[m]
#                     Δϵ = ϵₙ - ϵₘ

#                     fₙ = fermi(ϵₙ, μ, β)
#                     fₘ = fermi(ϵₘ, μ, β)
#                     Δf = (fₘ - fₙ)

#                     coupling = adiabatic_derivative[I][n,m] * adiabatic_derivative[J][m,n]
#                     if deriv
#                         calc.friction[I,J] += -2π * coupling * gauss(Δϵ, σ) * ∂fermi(ϵₘ, μ, β)
#                     else
#                         calc.friction[I,J] += 2π * coupling * gauss(Δϵ, σ) * Δf / Δϵ
#                     end
#                 end
#             end
#         end
#     end
# end
