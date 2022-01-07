
using NonadiabaticModels.DiabaticModels: LargeDiabaticModel

struct LargeDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::DependentField{Hermitian{T,Matrix{T}},Matrix{T}}
    derivative::DependentField{Matrix{Hermitian{T,Matrix{T}}},Matrix{T}}
    eigen::DependentField{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}},Matrix{T}}
    adiabatic_derivative::DependentField{Matrix{Matrix{T}},Matrix{T}}
    nonadiabatic_coupling::DependentField{Matrix{Matrix{T}},Matrix{T}}
    tmp_mat::Matrix{T}
    function LargeDiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        potential = Hermitian(zeros(n, n))
        derivative = [Hermitian(zeros(n, n)) for i=1:ndofs(model), j=1:atoms]
        eigen = Eigen(zeros(n), zeros(n,n)+I)
        adiabatic_derivative = [zeros(n, n) for i=1:ndofs(model), j=1:atoms]
        nonadiabatic_coupling = [zeros(n, n) for i=1:ndofs(model), j=1:atoms]
        tmp_mat = zeros(T, n, n)

        position = fill(NaN, ndofs(model), atoms)

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),
            tmp_mat,
        )
    end
end

function Calculator(model::LargeDiabaticModel, atoms::Integer, T::Type=Float64)
    LargeDiabaticCalculator{T}(model, atoms)
end

function evaluate_potential!(calc::LargeDiabaticCalculator, r)
    NonadiabaticModels.potential!(calc.model, calc.potential, r)
end

function evaluate_eigen!(calc::LargeDiabaticCalculator, r)
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
    eigen = get_eigen(calc, r)
    derivative = get_derivative(calc, r)
    for I in eachindex(derivative)
        LinearAlgebra.mul!(calc.tmp_mat, derivative[I], eigen.vectors)
        LinearAlgebra.mul!(calc.adiabatic_derivative[I], eigen.vectors', calc.tmp_mat)
    end
end

function evaluate_nonadiabatic_coupling!(calc::LargeDiabaticCalculator, r)
    eigen = get_eigen(calc, r)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)
    for I in eachindex(adiabatic_derivative)
        evaluate_nonadiabatic_coupling!(calc.nonadiabatic_coupling[I], adiabatic_derivative[I], eigen.values)
    end
end

function evaluate_nonadiabatic_coupling!(coupling::Matrix, adiabatic_derivative::Matrix, eigenvalues::Vector)
    for i=1:length(eigenvalues)
        for j=i+1:length(eigenvalues)
            coupling[j,i] = -adiabatic_derivative[j,i] / (eigenvalues[j]-eigenvalues[i])
            coupling[i,j] = -coupling[j,i]
        end
    end
end
