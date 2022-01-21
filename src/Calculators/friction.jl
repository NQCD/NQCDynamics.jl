
mutable struct FrictionCalculator{T,M} <: AbstractFrictionCalculator{T,M}
    model::M
    potential::T
    derivative::Matrix{T}
    friction::Matrix{T}
    function FrictionCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        new{T,M}(model, 0, zeros(ndofs(model), atoms), zeros(ndofs(model)*atoms, ndofs(model)*atoms))
    end
end

struct RingPolymerFrictionCalculator{T,M} <: AbstractFrictionCalculator{T,M}
    model::M
    potential::Vector{T}
    derivative::Array{T,3}
    friction::Array{T,3}
    function RingPolymerFrictionCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        new{T,M}(model, zeros(beads), zeros(ndofs(model), atoms, beads), zeros(ndofs(model)*atoms, ndofs(model)*atoms, beads))
    end
end

struct DiabaticFrictionCalculator{T,M} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::Hermitian{T,Matrix{T}}
    derivative::Matrix{Hermitian{T,Matrix{T}}}
    eigen::LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}}
    adiabatic_derivative::Matrix{Matrix{T}}
    nonadiabatic_coupling::Matrix{Matrix{T}}
    friction::Matrix{T}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function DiabaticFrictionCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        potential = Hermitian(zeros(n, n))
        derivative = [Hermitian(zeros(n, n)) for i=1:ndofs(model), j=1:atoms]
        eigen = Eigen(zeros(n), zeros(n, n)+I)
        adiabatic_derivative = [zeros(n, n) for i=1:ndofs(model), j=1:atoms]
        nonadiabatic_coupling = [zeros(n, n) for i=1:ndofs(model), j=1:atoms]
        friction = zeros(ndofs(model)*atoms, ndofs(model)*atoms)
        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T},n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)
        new{T,M}(model, potential, derivative, eigen,
                 adiabatic_derivative, nonadiabatic_coupling, friction,
                 tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
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
    NonadiabaticModels.friction!(calc.model, calc.friction, R)
end

function evaluate_friction!(calc::AbstractFrictionCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NonadiabaticModels.friction!(calc.model, calc.friction[:,:,i], R[:,:,i])
    end
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

    dofs = ndofs(calc.model)
    calc.friction .= 0
    for i in axes(R, 2) # Atoms
        for j in axes(R, 1) # dofs
            for m=2:nstates(calc.model)
                ω = calc.eigen.values[m] - calc.eigen.values[1]
                g = gauss(ω, calc.model.σ) / ω
                calc.friction[j+(i-1)*dofs] += 2π*abs2(calc.adiabatic_derivative[j,i][m,1])*g
            end
        end
    end
end
