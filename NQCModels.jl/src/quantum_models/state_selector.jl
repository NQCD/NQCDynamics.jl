using LinearAlgebra: eigen, Diagonal, diag

struct StateSelector{M<:QuantumModels.QuantumModel,B<:NQCBase.StateType} <: NQCModels.ClassicalModels.ClassicalModel
    quantum_model::M
    state::Int

    function StateSelector(quantum_model::QuantumModels.QuantumModel, state::Int, ::Type{B}) where {B<:NQCBase.StateType}
        state < 1 && throw(DomainError(state, "selected state must be greater than 0"))
        state > NQCModels.nstates(quantum_model) && throw(
            DomainError(state, "selected state must be less than or equal to the total number of states of the model"),
        )
        return new{typeof(quantum_model),B}(quantum_model, state)
    end
end

NQCModels.ndofs(model::StateSelector) = NQCModels.ndofs(model.quantum_model)

# generic out-of-place derivative, built on the in-place primitive below
function NQCModels.derivative(model::StateSelector, r::AbstractMatrix)
    output = zeros(eltype(r), size(r))
    NQCModels.derivative!(model, output, r)
    return output
end

# generic in-place potential — assumes a 1-element mutable container; see caveat above
function NQCModels.potential!(model::StateSelector, V::AbstractMatrix, r::AbstractMatrix)
    V .= NQCModels.potential(model, r)
    return V
end

# --- Diabatic model, want the adiabatic-basis energy of `state` ---
function NQCModels.potential(model::StateSelector{M,Adiabatic}, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    V = NQCModels.potential(model.quantum_model, r)
    return eigen(V).values[model.state]
end

function NQCModels.derivative!(model::StateSelector{M,Adiabatic}, output::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    V = NQCModels.potential(model.quantum_model, r)
    U = eigen(V).vectors
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = (U'*D[I]*U)[model.state, model.state]
    end
    return output
end

# --- Diabatic model, want the raw diagonal diabatic element ---
function NQCModels.potential(model::StateSelector{M,Diabatic}, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    V = NQCModels.potential(model.quantum_model, r)
    return V[model.state, model.state]
end

function NQCModels.derivative!(model::StateSelector{M,Diabatic}, output::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = D[I][model.state, model.state]
    end
    return output
end

# --- Adiabatic model, want the (already adiabatic) energy of `state` ---
function NQCModels.potential(model::StateSelector{M,Adiabatic}, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Adiabatic}}
    V = NQCModels.potential(model.quantum_model, r)
    return V[model.state]
end

function NQCModels.derivative!(model::StateSelector{M,Adiabatic}, output::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Adiabatic}}
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = D[I][model.state]
    end
    return output
end

struct ReducedQuantumModel{M<:QuantumModels.QuantumModel,B<:NQCBase.StateType} <: QuantumModels.QuantumModel{B}
    quantum_model::M
    states::Vector{Int}

    function ReducedQuantumModel(quantum_model::QuantumModels.QuantumModel, states::AbstractVector{Int}, ::Type{B}) where {B<:NQCBase.StateType}
        n = NQCModels.nstates(quantum_model)
        isempty(states) && throw(ArgumentError("`states` must not be empty"))
        any(s -> s < 1, states) && throw(DomainError(states, "selected states must be greater than 0"))
        any(s -> s > n, states) && throw(DomainError(states, "selected states must not exceed the total number of states ($n) of the underlying model"))
        allunique(states) || throw(ArgumentError("`states` must not contain duplicates"))
        return new{typeof(quantum_model),B}(quantum_model, collect(states))
    end
end

NQCModels.ndofs(model::ReducedQuantumModel) = NQCModels.ndofs(model.quantum_model)
NQCModels.nstates(model::ReducedQuantumModel) = length(model.states)

# generic out-of-place wrappers, built on the in-place primitives below
function NQCModels.potential(model::ReducedQuantumModel, r::AbstractMatrix)
    n = NQCModels.nstates(model)
    V = zeros(eltype(r), n, n)
    NQCModels.potential!(model, V, r)
    return V
end

function NQCModels.derivative(model::ReducedQuantumModel, r::AbstractMatrix)
    D_full = NQCModels.derivative(model.quantum_model, r)  # only used to get dof-array shape/eltype structure
    n = NQCModels.nstates(model)
    output = similar(D_full, Matrix{eltype(r)})
    for I in eachindex(output)
        output[I] = zeros(eltype(r), n, n)
    end
    NQCModels.derivative!(model, output, r)
    return output
end

# --- Diabatic underlying model, stay in diabatic representation ---
# Just take the submatrix — no diagonalization needed.
function NQCModels.potential!(model::ReducedQuantumModel{M,Diabatic}, V::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    Vfull = NQCModels.potential(model.quantum_model, r)
    V .= Vfull[model.states, model.states]
    return V
end

function NQCModels.derivative!(model::ReducedQuantumModel{M,Diabatic}, output::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = D[I][model.states, model.states]
    end
    return output
end

# --- Diabatic underlying model, truncate to a subset of adiabatic eigenstates ---
# Diagonalize the FULL Hamiltonian (truncating states before diagonalizing would
# corrupt the eigenbasis), then slice the requested rows/cols out of the result.
function NQCModels.potential!(model::ReducedQuantumModel{M,Adiabatic}, V::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    Vfull = NQCModels.potential(model.quantum_model, r)
    eigenvalues = eigen(Vfull).values
    V .= Diagonal(eigenvalues[model.states])
    return V
end

function NQCModels.derivative!(model::ReducedQuantumModel{M,Adiabatic}, output::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Diabatic}}
    V = NQCModels.potential(model.quantum_model, r)
    U = eigen(V).vectors
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = (U' * D[I] * U)[model.states, model.states]
    end
    return output
end

# --- Adiabatic underlying model, stay in adiabatic representation ---
function NQCModels.potential!(model::ReducedQuantumModel{M,Adiabatic}, V::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Adiabatic}}
    Vfull = NQCModels.potential(model.quantum_model, r)
    V .= Diagonal(Vfull[model.states])
    return V
end

function NQCModels.derivative!(model::ReducedQuantumModel{M,Adiabatic}, output::AbstractMatrix, r::AbstractMatrix) where {M<:QuantumModels.QuantumModel{Adiabatic}}
    D = NQCModels.derivative(model.quantum_model, r)
    for I in eachindex(output, D)
        output[I] = D[I][model.states, model.states]
    end
    return output
end