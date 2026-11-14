
struct AveragedPotential{M,T} <: ClassicalModel
    models::M
    tmp_derivative::Matrix{T}
end

AveragedPotential(models, r) = AveragedPotential(models, zero(r))

NQCModels.ndofs(model::AveragedPotential) = NQCModels.ndofs(model.models[1])

function NQCModels.potential(model::AveragedPotential, r::AbstractMatrix)
    V = sum(NQCModels.potential.(model.models, Ref(r)))
    return V / length(model.models)
end

function NQCModels.potential!(model::AveragedPotential, V::Matrix{<:Number}, r::AbstractMatrix)
    V .= sum(NQCModels.potential.(model.models, Ref(r))) / length(model.models)
end

function NQCModels.derivative(model::AveragedPotential, r::AbstractMatrix)
    D = zeros(size(model.tmp_derivative))
    for m in model.models
        NQCModels.derivative!(m, model.tmp_derivative, r)
        D .+= model.tmp_derivative
    end
    return D ./= length(model.models)
end

function NQCModels.derivative!(model::AveragedPotential, D::AbstractMatrix, r::AbstractMatrix)
    fill!(D, zero(eltype(D)))
    for m in model.models
        NQCModels.derivative!(m, model.tmp_derivative, r)
        D .+= model.tmp_derivative
    end
    D ./= length(model.models)
end
