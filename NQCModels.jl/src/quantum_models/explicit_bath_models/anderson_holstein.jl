
"""
## Newns-Anderson model
The **Anderson Impurity Model (AIM)** describes a localized impurity state interacting with a continuous band of bath states. AIM is a foundamental model in condensed matter physics and quantum chemistry, introduced by P.W. Anderson in 1961. The **Newns-Anderson model** is a generalization of the AIM, which includes the possibility of **multiple impurity states** and a more complex interaction with the bath. A key advantage of using the AIM lies in its ability to yield **analytical solutions** for the energy level distribution and the hybridization (coupling) density, making it a powerful tool for theoretical analysis.
"""
struct AndersonHolstein{M<:QuantumModel,B,T} <: QuantumModel
    impurity_model::M
    bath_model::B
    fermi_level::T
    nelectrons::Int
    couplings_rescale::Float64
    impurity_potential::Hermitian{Float64, Matrix{Float64}}
    impurity_derivative::Matrix{Hermitian{Float64, Matrix{Float64}}}
end

# include if statement here so that impurity_derivative always has shape Matrix{<:Hermitian} 
# that way it can be populated in derivative!(AndersonHolstien) and the shape is certain
function AndersonHolstein(impurity_model, bath; fermi_level=0.0, couplings_rescale=1.0) 
    fermi_level = austrip(fermi_level)
    nelectrons = count(bath.bathstates .≤ fermi_level)
    imp_potential = Hermitian(zeros(nstates(impurity_model),nstates(impurity_model)))
    imp_derivative = NQCModels.zero_derivative(impurity_model, hcat([0.0 for _ in NQCModels.dofs(impurity_model)]))
    return AndersonHolstein(impurity_model, bath, fermi_level, nelectrons, couplings_rescale, imp_potential, imp_derivative)
end

NQCModels.nstates(model::AndersonHolstein) = NQCModels.nstates(model.bath_model) + 1
NQCModels.ndofs(model::AndersonHolstein) = NQCModels.ndofs(model.impurity_model)
NQCModels.nelectrons(model::AndersonHolstein) = model.nelectrons
NQCModels.fermilevel(model::AndersonHolstein) = model.fermi_level

function NQCModels.potential!(model::AndersonHolstein, V::Hermitian, r::AbstractMatrix)
    NQCModels.potential!(model.impurity_model, model.impurity_potential, r)
    V[1,1] = model.impurity_potential[2,2] - model.impurity_potential[1,1]
    fillbathstates!(V, model.bath_model)
    fillbathcoupling!(V, model.impurity_potential[2,1], model.bath_model, model.couplings_rescale)

    return nothing
end

function NQCModels.derivative(model::AndersonHolstein, R::AbstractMatrix)
    D = [Hermitian(zero(matrix_template(model, eltype(R)))) for _=1:size(R, 1), _=1:size(R, 2)]
    @info "AndersonHolstein generated derivative size:" size = size(D)
    NQCModels.derivative!(model, D, R)
    return D
end

#= These multiple dispatch versions exits because the derivative of the 
impurity model may either be of type Hermitian or Matrix{Hermitian} 
depending on the number of spatial dimensions the impurity is defined over =#

"""
    function NQCModels.derivative!(model::AndersonHolstein, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
        output: nothing

Updates the derivative of the Anderson-Holstein Hamiltonian with respect to all spatial degrees of freedom to have the correct values
for a given position R.
    
This function is multiple dispatched over the shape of derivative(model.impurity_model) as these impurities may be defined over different
numbers of spatial degrees of freedom.
"""
function NQCModels.derivative!(model::AndersonHolstein, D::AbstractMatrix{<:Hermitian}, r::AbstractMatrix)
    # Get impurity model derivative
    NQCModels.derivative!(model.impurity_model, model.impurity_derivative, r)

    # Write and apply bath coupling
    @inbounds for i in axes(r, 1) #All model degrees of freedom
        @inbounds for j in axes(r, 2) # All particles
            D[i,j][1,1] = model.impurity_derivative[i,j][2,2] - model.impurity_derivative[i,j][1,1]
            fillbathcoupling!(D[i,j], model.impurity_derivative[i,j][2,1], model.bath_model, model.couplings_rescale)
        end
    end
    
    return nothing
end

function NQCModels.state_independent_potential(model::AndersonHolstein, r::AbstractMatrix)
    Vsystem = NQCModels.potential(model.impurity_model, r)
    return Vsystem[1,1]
end

function NQCModels.state_independent_potential!(model::AndersonHolstein, Vsystem::AbstractMatrix, r::AbstractMatrix)
    NQCModels.potential!(model.impurity_model, model.impurity_potential, r)
    @. Vsystem = model.impurity_potential
    return nothing
end

function NQCModels.state_independent_derivative(model::AndersonHolstein, r::AbstractMatrix)
    D_impurity_model = NQCModels.zero_derivative(model.impurity_model, r)
    NQCModels.derivative!(model.impurity_model, D_impurity_model, r)

    ∂V = zeros(size(r))
    @assert ∂V |> length == D_impurity_model |> length

    for I in eachindex(∂V, D_impurity_model)
        ∂V[I] = D_impurity_model[I][1,1]
    end

    return ∂V
end

function NQCModels.state_independent_derivative!(model::AndersonHolstein, ∂V::AbstractMatrix, r::AbstractMatrix)
    NQCModels.derivative!(model.impurity_model, model.impurity_derivative, r)
    
    for I in eachindex(∂V)
        ∂V[I] = model.impurity_derivative[I][1,1]
    end
end
