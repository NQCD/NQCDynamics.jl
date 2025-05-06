using NQCModels.DiabaticModels: LargeDiabaticModel

struct DensityMatrixCalculator{T,M} <: AbstractLargeDiabaticCalculators{T,M}
    model::M
    potential::DependentField{Hermitian{T,Matrix{T}},Matrix{T}}
    derivative::DependentField{Matrix{Hermitian{T,Matrix{T}}},Matrix{T}}
    eigen::DependentField{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}},Matrix{T}}
    adiabatic_derivative::DependentField{Matrix{Matrix{T}},Matrix{T}}
    nonadiabatic_coupling::DependentField{Matrix{Matrix{T}},Matrix{T}}
    tmp_mat::Matrix{T}
    stats::Dict{Symbol,Int}
end

function DensityMatrixCalculator{T}(model::M, atoms::Integer) where {T,M<:Model}
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

    return DensityMatrixCalculator{T,M}(
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



# ----------------------------------- from `large_diabatic.jl` ----------------------------------- #

# `get_adiabatic_derivative()` calls this function when updating
function evaluate_adiabatic_derivative!(calc::DensityMatrixCalculator, r)
    calc.stats[:adiabatic_derivative] += 1
    eigen = get_eigen(calc, r)
    derivative = get_derivative(calc, r) # calculate derivated of chosen model
    @inbounds for i in NQCModels.(calc)
        for j in NQCModels.dofs(calc)
            LinearAlgebra.mul!(calc.tmp_mat, derivative[j,i], eigen.vectors)
            LinearAlgebra.mul!(calc.adiabatic_derivative[j,i], eigen.vectors', calc.tmp_mat)
        end
    end
end

# `get_nonadiabatic_coupling()` calls this function when updating
function evaluate_nonadiabatic_coupling!(calc::DensityMatrixCalculator, r)
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
# ------------------------------------------------------------------------------------------------ #