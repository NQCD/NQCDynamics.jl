
using NonadiabaticModels.DiabaticModels: LargeDiabaticModel

struct LargeDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{T,M}
    model::M
    potential::Hermitian{T,Matrix{T}}
    derivative::Matrix{Hermitian{T,Matrix{T}}}
    eigen::LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}}
    adiabatic_derivative::Matrix{Matrix{T}}
    nonadiabatic_coupling::Matrix{Matrix{T}}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function LargeDiabaticCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
        n = nstates(model)
        potential = Hermitian(zeros(n, n))
        derivative = [Hermitian(zeros(n, n)) for i=1:ndofs(model), j=1:atoms]
        eigen = Eigen(zeros(n), zeros(n,n)+I)
        adiabatic_derivative = [zeros(n, n) for i=1:ndofs(model), j=1:atoms]
        nonadiabatic_coupling = [zeros(n, n) for i=1:ndofs(model), j=1:atoms]
        tmp_mat = zeros(T, n, n)
        tmp_mat_complex1 = zeros(Complex{T}, n, n)
        tmp_mat_complex2 = zeros(Complex{T}, n, n)
        new{T,M}(model, potential, derivative, eigen, adiabatic_derivative, nonadiabatic_coupling,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

function Calculator(model::LargeDiabaticModel, atoms::Integer, T::Type=Float64)
    LargeDiabaticCalculator{T}(model, atoms)
end

function evaluate_potential!(calc::LargeDiabaticCalculator, R)
    NonadiabaticModels.potential!(calc.model, calc.potential, R)
end

function evaluate_potential!(calc::LargeDiabaticCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NonadiabaticModels.potential!(calc.model, calc.potential[i], R[:,:,i])
    end
end

function transform_derivative!(calc::LargeDiabaticCalculator)
    for I in eachindex(calc.derivative)
        LinearAlgebra.mul!(calc.tmp_mat, calc.derivative[I], calc.eigen.vectors)
        LinearAlgebra.mul!(calc.adiabatic_derivative[I], calc.eigen.vectors', calc.tmp_mat)
    end
end

function eigen!(calc::LargeDiabaticCalculator)
    eig = LinearAlgebra.eigen(calc.potential)
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

function evaluate_nonadiabatic_coupling!(calc::LargeDiabaticCalculator)
    for I in eachindex(calc.adiabatic_derivative)
        evaluate_nonadiabatic_coupling!(calc.nonadiabatic_coupling[I], calc.adiabatic_derivative[I], calc.eigen.values)
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
