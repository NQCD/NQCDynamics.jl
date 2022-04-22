
struct RingPolymerLargeDiabaticCalculator{T,M} <: AbstractLargeDiabaticCalculators{T,M}
    model::M
    potential::DependentField{Vector{Hermitian{T,Matrix{T}}},Array{T,3}}
    derivative::DependentField{Array{Hermitian{T,Matrix{T}},3},Array{T,3}}
    eigen::DependentField{Vector{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}}},Array{T,3}}
    adiabatic_derivative::DependentField{Array{Matrix{T},3},Array{T,3}}
    nonadiabatic_coupling::DependentField{Array{Matrix{T},3},Array{T,3}}
    tmp_mat::Matrix{T}
    stats::Dict{Symbol,Int}
    function RingPolymerLargeDiabaticCalculator{T}(model::M, atoms::Integer, beads::Integer) where {T,M<:Model}
        mat = NQCModels.DiabaticModels.matrix_template(model, T)
        vec = NQCModels.DiabaticModels.vector_template(model, T)

        potential = [Hermitian(zero(mat)) for _=1:beads]
        derivative = [Hermitian(zero(mat)) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        eigen = [Eigen(zero(vec), zero(mat)+I) for _=1:beads]
        adiabatic_derivative = [zero(mat) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        nonadiabatic_coupling = [zero(mat) for _=1:ndofs(model), _=1:atoms, _=1:beads]
        tmp_mat = zero(mat)

        position = fill(NaN, ndofs(model), atoms, beads)

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

function Calculator(model::LargeDiabaticModel, atoms::Integer, beads::Integer, T::Type=Float64)
    RingPolymerLargeDiabaticCalculator{T}(model, atoms, beads)
end

function evaluate_potential!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:potential] += 1
    @views @inbounds for i in beads(calc)
        NQCModels.potential!(calc.model, calc.potential[i], r[:,:,i])
    end
end

function evaluate_eigen!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:eigen] += 1
    potential = get_potential(calc, r)
    @inbounds for i in beads(calc)
        eig = LinearAlgebra.eigen(potential[i])
        correct_phase!(eig, calc.eigen[i].vectors)
        copy!(calc.eigen[i].vectors, eig.vectors)
        copy!(calc.eigen[i].values, eig.values)
    end
end

function evaluate_adiabatic_derivative!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:adiabatic_derivative] += 1
    eigen = get_eigen(calc, r)
    derivative = get_derivative(calc, r)
    @inbounds for i in beads(calc)
        for j in mobileatoms(calc)
            for k in dofs(calc)
                LinearAlgebra.mul!(calc.tmp_mat, derivative[k,j,i], eigen[i].vectors)
                LinearAlgebra.mul!(calc.adiabatic_derivative[k,j,i], eigen[i].vectors', calc.tmp_mat)
            end
        end
    end
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:nonadiabatic_coupling] += 1
    eigen = get_eigen(calc, r)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)

    @inbounds for i in beads(calc)
        evaluate_inverse_difference_matrix!(calc.tmp_mat, eigen[i].values)
        for j in mobileatoms(calc)
            for k in dofs(calc)
                multiply_elementwise!(calc.nonadiabatic_coupling[k,j,i], adiabatic_derivative[k,j,i], calc.tmp_mat)
            end
        end
    end
end
