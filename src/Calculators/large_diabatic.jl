
struct LargeDiabaticCalculator{T,M} <: AbstractDiabaticCalculator{M}
    model::M
    potential::Hermitian{T,Matrix{T}}
    derivative::Matrix{Hermitian{T,Matrix{T}}}
    eigenvalues::Vector{T}
    eigenvectors::Matrix{T}
    adiabatic_derivative::Matrix{Matrix{T}}
    nonadiabatic_coupling::Matrix{Matrix{T}}
    tmp_mat::Matrix{T}
    tmp_mat_complex1::Matrix{Complex{T}}
    tmp_mat_complex2::Matrix{Complex{T}}
    function LargeDiabaticCalculator{T}(model::M, DoFs::Integer, atoms::Integer) where {T,M<:Model}
        potential = Hermitian(zeros(model.n_states, model.n_states))
        derivative = [Hermitian(zeros(model.n_states, model.n_states)) for i=1:DoFs, j=1:atoms]
        eigenvalues = zeros(model.n_states)
        eigenvectors = zeros(model.n_states, model.n_states)
        adiabatic_derivative = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        nonadiabatic_coupling = [zeros(model.n_states, model.n_states) for i=1:DoFs, j=1:atoms]
        tmp_mat = zeros(T, model.n_states, model.n_states)
        tmp_mat_complex1 = zeros(Complex{T}, model.n_states, model.n_states)
        tmp_mat_complex2 = zeros(Complex{T}, model.n_states, model.n_states)
        new{T,M}(model, potential, derivative, eigenvalues, eigenvectors, adiabatic_derivative, nonadiabatic_coupling,
            tmp_mat, tmp_mat_complex1, tmp_mat_complex2)
    end
end

function Calculator(model::LargeDiabaticModel, DoFs::Integer, atoms::Integer, T::Type=Float64)
    LargeDiabaticCalculator{T}(model, DoFs, atoms)
end

evaluate_potential!(calc::LargeDiabaticCalculator, R) = potential!(calc.model, calc.potential, R)

function evaluate_potential!(calc::LargeDiabaticCalculator, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        potential!(calc.model, calc.potential[i], R[:,:,i])
    end
end

function transform_derivative!(calc::LargeDiabaticCalculator)
    for I in eachindex(calc.derivative)
        mul!(calc.tmp_mat, calc.derivative[I], calc.eigenvectors)
        mul!(calc.adiabatic_derivative[I], calc.eigenvectors', calc.tmp_mat)
    end
end

function eigen!(calc::LargeDiabaticCalculator)
    eig = eigen(calc.potential)
    correct_phase!(eig, calc.eigenvectors)
    copyto!(calc.eigenvectors, eig.vectors)
    copyto!(calc.eigenvalues, eig.values)
end

function correct_phase!(eig::Eigen, old_eigenvectors::AbstractMatrix)
    @views for i=1:length(eig.values)
        if dot(eig.vectors[:,i], old_eigenvectors[:,i]) < 0
            eig.vectors[:,i] .*= -1
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::LargeDiabaticCalculator)
    for I in eachindex(calc.adiabatic_derivative)
        evaluate_nonadiabatic_coupling!(calc.nonadiabatic_coupling[I], calc.adiabatic_derivative[I], calc.eigenvalues)
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
