
struct RingPolymerLargeDiabaticCalculator{T,M} <: AbstractLargeDiabaticCalculators{T,M}
    model::M
    potential::DependentField{Vector{Hermitian{T,Matrix{T}}},Array{T,3}}
    derivative::DependentField{Array{Hermitian{T,Matrix{T}},3},Array{T,3}}
    eigen::DependentField{Vector{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}}},Array{T,3}}
    adiabatic_derivative::DependentField{Array{Matrix{T},3},Array{T,3}}
    nonadiabatic_coupling::DependentField{Array{Matrix{T},3},Array{T,3}}
    tmp_mat::Matrix{T}

    centroid::DependentField{Matrix{T},Array{T,3}}
    centroid_potential::DependentField{Hermitian{T,Matrix{T}},Array{T,3}}
    centroid_derivative::DependentField{Matrix{Hermitian{T,Matrix{T}}},Array{T,3}}
    centroid_eigen::DependentField{LinearAlgebra.Eigen{T,T,Matrix{T},Vector{T}},Array{T,3}}
    centroid_adiabatic_derivative::DependentField{Matrix{Matrix{T}},Array{T,3}}
    centroid_nonadiabatic_coupling::DependentField{Matrix{Matrix{T}},Array{T,3}}

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

        centroid = zeros(T, ndofs(model), atoms)
        centroid_potential = Hermitian(zero(mat))
        centroid_derivative = [Hermitian(zero(mat)) for _=1:ndofs(model), _=1:atoms]
        centroid_eigen = Eigen(zero(vec), zero(mat) + I)
        centroid_adiabatic_derivative = [zero(mat) for _ in CartesianIndices(centroid_derivative)]
        centroid_nonadiabatic_coupling = [zero(mat) for _ in CartesianIndices(centroid_derivative)]

        position = fill(NaN, ndofs(model), atoms, beads)

        stats = Dict{Symbol,Int}(
            :potential=>0,
            :derivative=>0,
            :eigen=>0,
            :adiabatic_derivative=>0,
            :nonadiabatic_coupling=>0,
            :centroid=>0,
            :centroid_potential=>0,
            :centroid_derivative=>0,
            :centroid_eigen=>0,
            :centroid_adiabatic_derivative=>0,
            :centroid_nonadiabatic_coupling=>0,
        )

        new{T,M}(
            model,
            DependentField(potential, copy(position)),
            DependentField(derivative, copy(position)),
            DependentField(eigen, copy(position)),
            DependentField(adiabatic_derivative, copy(position)),
            DependentField(nonadiabatic_coupling, copy(position)),
            tmp_mat,

            DependentField(centroid, copy(position)),
            DependentField(centroid_potential, copy(position)),
            DependentField(centroid_derivative, copy(position)),
            DependentField(centroid_eigen, copy(position)),
            DependentField(centroid_adiabatic_derivative, copy(position)),
            DependentField(centroid_nonadiabatic_coupling, copy(position)),

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

function evaluate_centroid_potential!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:centroid_potential] += 1
    centroid = get_centroid(calc, r)
    NQCModels.potential!(calc.model, calc.centroid_potential, centroid)
    return nothing
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

function evaluate_centroid_eigen!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:centroid_eigen] += 1
    potential = get_centroid_potential(calc, r)
    eig = LinearAlgebra.eigen(potential)
    correct_phase!(eig, calc.centroid_eigen.vectors)
    copy!(calc.centroid_eigen.vectors, eig.vectors)
    copy!(calc.centroid_eigen.values, eig.values)
    return nothing
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

function evaluate_centroid_adiabatic_derivative!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:centroid_adiabatic_derivative] += 1
    eigen = get_centroid_eigen(calc, r)
    derivative = get_centroid_derivative(calc, r)
    @inbounds for j in mobileatoms(calc)
        for k in dofs(calc)
            LinearAlgebra.mul!(calc.tmp_mat, derivative[k,j], eigen.vectors)
            LinearAlgebra.mul!(calc.centroid_adiabatic_derivative[k,j], eigen.vectors', calc.tmp_mat)
        end
    end
    return nothing
end

function evaluate_nonadiabatic_coupling!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:nonadiabatic_coupling] += 1
    eigen = get_eigen(calc, r)
    adiabatic_derivative = get_adiabatic_derivative(calc, r)

    @inbounds for i in beads(calc)
        evaluate_inverse_difference_matrix!(calc.tmp_mat, eigen[i].values)
        for j in mobileatoms(calc)
            for k in dofs(calc)
                @. calc.nonadiabatic_coupling[k,j,i] = adiabatic_derivative[k,j,i] * calc.tmp_mat
            end
        end
    end
end

function evaluate_centroid_nonadiabatic_coupling!(calc::RingPolymerLargeDiabaticCalculator, r::AbstractArray{T,3}) where {T}
    calc.stats[:centroid_nonadiabatic_coupling] += 1
    adiabatic_derivative = get_centroid_adiabatic_derivative(calc, r)
    eigen = get_centroid_eigen(calc, r)
    evaluate_inverse_difference_matrix!(calc.tmp_mat, eigen.values)
    @inbounds for j in mobileatoms(calc)
        for k in dofs(calc)
            @. calc.centroid_nonadiabatic_coupling[j,k] = adiabatic_derivative[j,k] * calc.tmp_mat
        end
    end
end
