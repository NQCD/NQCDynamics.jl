
"""
	TensorialFriction
	
Abstract model type which yields full-rank friction matrices. 

If atomic masses are required to calculate friction in your ElectronicFrictionProvider (e.g. for Isotope support), the atomic masses to use should be included as a type field. 

A `friction_atoms` field must be included - This indicates that the friction matrix is only partially returned. 
Subtypes of this must implement `NQCModels.FrictionModels.get_friction_matrix(model::ACEdsODF, R::AbstractMatrix, friction_atoms::AbstractVector) --> ::AbstractMatrix{eltype(R)}`.
Subtypes of this must implement `ndofs(m::TensorialFriction) --> ::Int`.

Units of friction are mass-weighted, and the atomic unit of friction is: [E_h / ħ / m_e]
Convert common friction units such as ps^-1 or meV ps Å^-2 using `UnitfulAtomic.austrip`. 

A system-size friction matrix is obtained with `friction(model, positions)`, or `friction!(model, friction_matrix, positions)`
"""
abstract type TensorialFriction <: ElectronicFrictionProvider end

function get_friction_matrix(model::ConstantFriction{AbstractMatrix}, ::AbstractMatrix)
	return model.γ
end

function get_friction_matrix(model::RandomFriction, R::AbstractMatrix)
	F = randn(model.ndofs * size(R, 2), model.ndofs * size(R, 2))
	F .= F'F
	F .= (F + F')/2
end

function friction(model::TensorialFriction, positions::AbstractMatrix)
	friction_matrix = similar(positions, [*(size(positions)...) for i in 1:2]...) # Generate a (ndofs * n_atoms)² matrix to hold friction
	friction!(model, friction_matrix, positions)
	return friction_matrix
end

function friction!(model::TensorialFriction, F::AbstractMatrix, R::AbstractMatrix)
	indices = friction_matrix_indices(model.friction_atoms, NQCModels.ndofs(model))
	F[indices, indices] .= get_friction_matrix(model, R)
end


