
using NQCModels.DiabaticModels: LargeDiabaticModel

abstract type AbstractLargeDiabaticCalculators{T,M} <: AbstractDiabaticCalculator{T,M} end

struct LargeDiabaticCalculator{T,M} <: AbstractLargeDiabaticCalculators{T,M}
    model::M
    potential::DependentField{Hermitian{T,Matrix{T}},Matrix{T}}
    derivative::DependentField{Matrix{Hermitian{T,Matrix{T}}},Matrix{T}}
    eigen::DependentField{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}},Matrix{T}}
    adiabatic_derivative::DependentField{Matrix{Matrix{T}},Matrix{T}}
    nonadiabatic_coupling::DependentField{Matrix{Matrix{T}},Matrix{T}}
    tmp_mat::Matrix{T}
    stats::Dict{Symbol,Int}
    function LargeDiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        mat = NQCModels.DiabaticModels.matrix_template(model, T)
        vec = NQCModels.DiabaticModels.vector_template(model, T)

        potential = Hermitian(zero(mat))
        derivative = [Hermitian(zero(mat)) for i=1:ndofs(model), j=1:atoms]
        eigen = Eigen(zero(vec), zero(mat)+I)
        adiabatic_derivative = [zero(mat) for i=1:ndofs(model), j=1:atoms]
        nonadiabatic_coupling = [zero(mat) for i=1:ndofs(model), j=1:atoms]
        tmp_mat = zero(mat)

        position = fill(NaN, ndofs(model), atoms)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
            :eigen=>0,
            :adiabatic_derivative=>0,
            :nonadiabatic_coupling=>0
        )

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),
            tmp_mat,
            stats,
        )
    end
end

function Calculator(model::LargeDiabaticModel, atoms::Integer, T::Type=Float64)
    LargeDiabaticCalculator{T}(model, atoms)
end

function evaluate_potential!(calc::LargeDiabaticCalculator, r)
    calc.stats[:potential] += 1
    NQCModels.potential!(calc.model, calc.potential, r)
    return nothing
end

function evaluate_eigen!(calc::LargeDiabaticCalculator, r)
    calc.stats[:eigen] += 1
    potential = get_potential(calc, r)
    eig = LinearAlgebra.eigen(potential)
    correct_phase!(eig, calc.eigen.vectors)
    copyto!(calc.eigen.vectors, eig.vectors)
    copyto!(calc.eigen.values, eig.values)
end

function correct_phase!(eig::LinearAlgebra.Eigen, old_eigenvectors::AbstractMatrix)
    @views for i=1:length(eig.values)
        if LinearAlgebra.dot(eig.vectors[:,i], old_eigenvectors[:,i]) < 0
            eig.vectors[:,i] .*= -1
        end
    end
end

function evaluate_adiabatic_derivative!(calc::LargeDiabaticCalculator, r)
    calc.stats[:adiabatic_derivative] += 1
    eigen = get_eigen(calc, r)
    derivative = get_derivative(calc, r)
    @inbounds for i in NQCModels.mobileatoms(calc)
        for j in NQCModels.dofs(calc)
            LinearAlgebra.mul!(calc.tmp_mat, derivative[j,i], eigen.vectors)
            LinearAlgebra.mul!(calc.adiabatic_derivative[j,i], eigen.vectors', calc.tmp_mat)
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::LargeDiabaticCalculator, r)
    calc.stats[:nonadiabatic_coupling] += 1
    eigen = get_eigen(calc, r)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)

    evaluate_inverse_difference_matrix!(calc.tmp_mat, eigen.values)

    @inbounds for i in NQCModels.mobileatoms(calc)
        for j in NQCModels.dofs(calc)
            multiply_elementwise!(calc.nonadiabatic_coupling[j,i], adiabatic_derivative[j,i], calc.tmp_mat)
        end
    end
end

function evaluate_inverse_difference_matrix!(out, eigenvalues)
    @inbounds for i in eachindex(eigenvalues)
        for j in eachindex(eigenvalues)
            out[j,i] = 1 / (eigenvalues[i] - eigenvalues[j])
        end
        out[i,i] = zero(eltype(out))
    end
end

function multiply_elementwise!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalue_difference_matrix::Matrix)
    @inbounds for I in eachindex(coupling, adiabatic_derivative, eigenvalue_difference_matrix)
        coupling[I] = adiabatic_derivative[I] * eigenvalue_difference_matrix[I]
    end
end
