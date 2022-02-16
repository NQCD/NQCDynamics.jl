
mutable struct FrictionCalculator{T,M} <: AbstractFrictionCalculator{T,M}
    model::M
    potential::DependentField{T,Matrix{T}}
    derivative::DependentField{Matrix{T},Matrix{T}}
    friction::DependentField{Matrix{T},Matrix{T}}
    function FrictionCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}

        potential = zero(T)
        derivative = zeros(ndofs(model), atoms)
        friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms)

        position = fill(NaN, ndofs(model), atoms)

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(friction, copy(position)),
        )
    end
end

struct RingPolymerFrictionCalculator{T,M} <: AbstractFrictionCalculator{T,M}
    model::M
    potential::DependentField{Vector{T},Array{T,3}}
    derivative::DependentField{Array{T,3},Array{T,3}}
    friction::DependentField{Array{T,3},Array{T,3}}
    function RingPolymerFrictionCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}

        potential = zeros(beads)
        derivative = zeros(ndofs(model), atoms, beads)
        friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms, beads)

        position = fill(NaN, ndofs(model), atoms, beads)

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(friction, copy(position)),
        )
    end
end

struct DiabaticFrictionCalculator{T,M} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::DependentField{Hermitian{T,Matrix{T}},Matrix{T}}
    derivative::DependentField{Matrix{Hermitian{T,Matrix{T}}},Matrix{T}}
    eigen::DependentField{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}},Matrix{T}}
    adiabatic_derivative::DependentField{Matrix{Matrix{T}},Matrix{T}}
    nonadiabatic_coupling::DependentField{Matrix{Matrix{T}},Matrix{T}}
    friction::DependentField{Matrix{T},Matrix{T}}
    tmp_mat::Matrix{T}

    function DiabaticFrictionCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        mat = NQCModels.DiabaticModels.matrix_template(model, T)
        vec = NQCModels.DiabaticModels.vector_template(model, T)

        potential = Hermitian(zero(mat))
        derivative = [Hermitian(zero(mat)) for i=1:ndofs(model), j=1:atoms]
        eigen = Eigen(zero(vec), zero(mat)+I)
        adiabatic_derivative = [zero(mat) for i=1:ndofs(model), j=1:atoms]
        nonadiabatic_coupling = [zero(mat) for i=1:ndofs(model), j=1:atoms]
        friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms)
        tmp_mat = zero(mat)

        position = fill(NaN, ndofs(model), atoms)

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),
            DependentField(friction, copy(position)),
            tmp_mat
        )
    end
end
function Calculator(model::AdiabaticFrictionModel, atoms::Integer, T::Type=Float64)
    FrictionCalculator{T}(model, atoms)
end
function Calculator(model::DiabaticFrictionModel, atoms::Integer, T::Type=Float64)
    DiabaticFrictionCalculator{T}(model, atoms)
end
function Calculator(model::AdiabaticFrictionModel, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerFrictionCalculator{T}(model, atoms, beads)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractMatrix)
    NQCModels.friction!(calc.model, calc.friction, R)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NQCModels.friction!(calc.model, calc.friction[:,:,i], R[:,:,i])
    end
end

@doc raw"""
    evaluate_friction!(calc::DiabaticFrictionCalculator, R::AbstractMatrix)

Evaluate the electronic friction for a model given in the diabatic representation.

```math
γ = 2πħ ∑ⱼ <1|dH|j><j|dH|1> δ(ωⱼ) / ωⱼ
```
Note that the delta function is approximated by a normalised gaussian.
"""
function evaluate_friction!(calc::DiabaticFrictionCalculator, r::AbstractMatrix)

    gauss(x, σ) = exp(-0.5 * x^2 / σ^2) / (σ*sqrt(2π))
    adiabatic_derivative = get_adiabatic_derivative(calc, r)
    eigen = get_eigen(calc, r)

    dofs = ndofs(calc.model)
    calc.friction .= 0
    for i in axes(r, 2) # Atoms
        for j in axes(r, 1) # dofs
            for m=2:nstates(calc.model)
                ω = eigen.values[m] - eigen.values[1]
                g = gauss(ω, calc.model.σ) / ω
                calc.friction[j+(i-1)*dofs] += 2π*abs2(adiabatic_derivative[j,i][m,1])*g
            end
        end
    end
end
