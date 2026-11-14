using LinearAlgebra

"""
    struct: WideBandBath <: QuantumFrictionModel

An explicit impurity-bath model which couples a quantum model to a discritised bath. The bath density of states
is explicitly given by the wide-band approximation and the bath is discretised using the trapezoidal rule.

It is recommended that this is the explicit bath model used in conjuction with MDEF and the friction evaluators defined in
NQCCalculators.jl. For the modelling of explicit bathstate dynamics (i.e. using dynamics methods such as IESH) it is
recommended that the user use the Anderson-Holstein model alongside the Shenvi-Gauss Legendre discretisation scheme.
"""
struct WideBandBath{M<:QuantumModel,V<:AbstractVector,T,D} <: QuantumFrictionModel
    impurity_model::M
    bathstates::V
    ρ::T
    tmp_derivative::Base.RefValue{D}
end

function WideBandBath(model::QuantumModel; step, bandmin, bandmax)
    bathstates = austrip.(collect(range(bandmin, bandmax; step=step)))
    ρ = length(bathstates) / (bathstates[end] - bathstates[begin])

    tmp_derivative = Ref(NQCModels.zero_derivative(model, zeros(1,1)))

    return WideBandBath(model, bathstates, ρ, tmp_derivative)
end

NQCModels.nstates(model::WideBandBath) = NQCModels.nstates(model.impurity_model) + length(model.bathstates) - 1
NQCModels.ndofs(model::WideBandBath) = NQCModels.ndofs(model.impurity_model)
NQCModels.nelectrons(model::WideBandBath) = fld(NQCModels.nstates(model), 2)
NQCModels.fermilevel(::WideBandBath) = 0.0

function NQCModels.potential(model::WideBandBath, r::AbstractMatrix)

    Vsystem = NQCModels.potential(model.impurity_model, r)
    n = NQCModels.nstates(model.impurity_model)
    ϵ0 = Vsystem[1,1]
    V = Hermitian(zeros((NQCModels.nstates(model), NQCModels.nstates(model))))

    # System states
    @views V[diagind(V)[begin:n-1]] .= Vsystem[diagind(Vsystem)[begin+1:end]] .- ϵ0

    # Bath states
    V[diagind(V)[n:end]] .= model.bathstates

    # Coupling
    for i=1:n-1
        V.data[i,n:end] .= Vsystem[i+1,1] / sqrt(model.ρ)
    end

    return V
end

function NQCModels.potential!(model::WideBandBath, V::Hermitian, r::AbstractMatrix)

    Vsystem = NQCModels.potential(model.impurity_model, r)
    n = NQCModels.nstates(model.impurity_model)
    ϵ0 = Vsystem[1,1]

    # System states
    @views V[diagind(V)[begin:n-1]] .= Vsystem[diagind(Vsystem)[begin+1:end]] .- ϵ0

    # Bath states
    V[diagind(V)[n:end]] .= model.bathstates

    # Coupling
    for i=1:n-1
        V.data[i,n:end] .= Vsystem[i+1,1] / sqrt(model.ρ)
    end
end


function NQCModels.derivative(model::WideBandBath, R::AbstractMatrix)
    D = [Hermitian(zero(matrix_template(model, eltype(R)))) for _=1:size(R, 1), _=1:size(R, 2)]
    NQCModels.derivative!(model, D, R)
    return D
end


"""
    function NQCModels.derivative!(model::WideBandBath, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
        output: nothing

Updates the derivative of the Anderson-Holstein Hamiltonian with respect to all spatial degrees of freedom to have the correct values
for a given position R.

This fucntion is multiple dispatched over the shape of derivative(model.impurity_model) as these sub-models may be defined over different
numbers of spatial degrees of freedom.
"""
function NQCModels.derivative!(model::WideBandBath, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)

    # Get inner model derivative
    Dsystem = NQCModels.zero_derivative(model.impurity_model, r)
    NQCModels.derivative!(model.impurity_model, Dsystem, r)
    n = NQCModels.nstates(model.impurity_model)

    for dof in axes(D, 1)
        for atom in axes(D, 2)
            ∂ϵ0 = Dsystem[dof, atom][1,1]

            # System states
            D_system_output = @view D[dof, atom][diagind(D[dof,atom])[begin:n-1]]
            D_system_input = @view Dsystem[dof, atom][diagind(Dsystem[dof, atom])[begin+1:end]]
            @. D_system_output = D_system_input - ∂ϵ0

            # Coupling
            for i=1:n-1
                D[dof,atom].data[i,n:end] .= Dsystem[dof,atom][i+1,1] / sqrt(model.ρ)
            end
        end
    end
end

function NQCModels.state_independent_potential(model::WideBandBath, r::AbstractMatrix)
    Vsystem = NQCModels.potential(model.impurity_model, r)
    return Vsystem[1,1]
end

function NQCModels.state_independent_potential!(model::WideBandBath, Vsystem::AbstractMatrix, r::AbstractMatrix)
    Vsystem .= NQCModels.potential(model.impurity_model, r)
end

function NQCModels.state_independent_derivative(model::WideBandBath, r::AbstractMatrix)
    ∂V = zeros(size(r))
    NQCModels.state_independent_derivative!(model, ∂V, r)
    return ∂V
end

function NQCModels.state_independent_derivative!(model::WideBandBath, ∂V::AbstractMatrix, r::AbstractMatrix)
    Dsystem = NQCModels.zero_derivative(model.impurity_model, r)
    NQCModels.derivative!(model.impurity_model, Dsystem, r)

    @assert ∂V |> length == Dsystem |> length

    for I in eachindex(∂V)
        ∂V[I] = Dsystem[I][1,1]
    end
end

function get_subsystem_derivative(model::WideBandBath, r::AbstractMatrix)
    if size(r) != size(model.tmp_derivative[])
        model.tmp_derivative[] = NQCModels.zero_derivative(model.impurity_model, r)
    end
    NQCModels.derivative!(model.impurity_model, model.tmp_derivative[], r)
    return model.tmp_derivative[]
end
