
struct AdiabaticStateSelector{M} <: NQCModels.ClassicalModels.ClassicalModel
    quantum_model::M
    state::Int
    function AdiabaticStateSelector(quantum_model, state)
        state < 1 && throw(DomainError(state, "selected state must be greater than 0"))
        state > NQCModels.nstates(quantum_model) && throw(
            DomainError(
                state,
                "selected state must be less than the total number of states of the diabatic model",
            ),
        )
        return new{typeof(quantum_model)}(quantum_model, state)
    end
end

NQCModels.ndofs(model::AdiabaticStateSelector) = NQCModels.ndofs(model.quantum_model)

function NQCModels.potential(model::AdiabaticStateSelector, r::AbstractMatrix)
    V = NQCModels.potential(model.quantum_model, r)
    eigenvalues = eigen(V).values
    return eigenvalues[model.state]
end

function NQCModels.potential!(model::AdiabaticStateSelector, V::Real, r::AbstractMatrix)
    Vsystem = NQCModels.potential(model.quantum_model, r)
    V .= eigen(Vsystem).values[model.state]
end

function NQCModels.derivative!(
    model::AdiabaticStateSelector,
    output::AbstractMatrix,
    r::AbstractMatrix,
)
    V = NQCModels.potential(model.quantum_model, r)
    U = eigen(V).vectors
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = (U'*D[I]*U)[model.state, model.state]
    end
    return output
end
