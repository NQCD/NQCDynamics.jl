
export Subsystem, CompositeModel
using .FrictionModels
using .ClassicalModels

"""
Subsystem(M, indices)

A subsystem is a Model which only applies to a subset of the degrees of freedom of the original model. 

**When combined in a CompositeModel**, `potential()`, `derivative!()` and `friction!()` will be sourced from the respective Subsystems. 

Calling `potential()`, `derivative!()`, or `friction!()` on a subsystem directly will output the respective values **for the entire system**. 

The Model specified will be supplied with the positions of the entire system for evaluation. 

"""
struct Subsystem{M<:Union{Model, FrictionModels.ElectronicFrictionProvider}}
	model::M
	indices::Union{Vector{Int}, Colon}
end

function Base.show(io::IO, subsystem::Subsystem)
    print(io, "Subsystem:\n\t🏎️ $(subsystem.model)\n\t🔢 $(subsystem.indices)\n")
end

function Subsystem(model, indices::Union{Int, UnitRange{Int}}=:)
	# Convert indices to a Vector{Int} or : for consistency
	if isa(indices, Int)
		indices = [indices:indices]
	elseif isa(indices, UnitRange{Int})
		indices = collect(indices)
	end
	Subsystem(model, indices)
end

# Passthrough functions to Model functions
potential(subsystem::Subsystem, R::AbstractMatrix) = potential(subsystem.model, R)
derivative(subsystem::Subsystem, R::AbstractMatrix) = derivative(subsystem.model, R)
derivative!(subsystem::Subsystem, D::AbstractMatrix, R::AbstractMatrix) = derivative!(subsystem.model, D, R)
FrictionModels.friction(subsystem::Subsystem, R::AbstractMatrix) = friction(subsystem.model, R)
FrictionModels.friction!(subsystem::Subsystem, F::AbstractMatrix, R::AbstractMatrix) = friction!(subsystem.model, F, R)
dofs(subsystem::Subsystem) = dofs(subsystem.model)
ndofs(subsystem::Subsystem) = ndofs(subsystem.model)

"""
CompositeModel(Subsystems...)

A CompositeModel is composed of multiple Subsystems, creating an effective model which evaluates each Subsystem for its respective indices. 
"""
struct CompositeModel{S<:Vector{<:Subsystem}, D<:Int} <: ClassicalModels.ClassicalModel
	subsystems::S
	ndofs::D
end

function Base.show(io::IO, model::CompositeModel)
    print(io, "CompositeModel with subsystems:\n", [system for system in model.subsystems]...)
end

"""
    CompositeModel(subsystems::Subsystem...)

Combine multiple Subsystems into a single model to be handled by NQCDynamics.jl in simulations. 
Any calls made to `potential`, `derivative` and `friction` will apply each subsystem's model to the respective atoms while ignoring any other atoms. 

Some checks are made to ensure each atom is affected by a model and that each model is applied over the same degrees of freedom, but no other sanity checks are made. 
"""
CompositeModel(subsystems::Subsystem...) = CompositeModel(check_models(subsystems...)...) # Check subsystems are a valid combination
# Convenience version which uses Models and automatically applies across all indices. 
CompositeModel(models::Model...) = CompositeModel(Subsystem.(models, Colon())...)
# Mixed use version with Models and Subsystems. 
function CompositeModel(model_or_subsystems::Union{Model, FrictionModels.ElectronicFrictionProvider, Subsystem}...)
	subsystems = Subsystem[]
	for model_or_subsystem in model_or_subsystems
		if !isa(model_or_subsystem, Subsystem)
			push!(subsystems, Subsystem(model_or_subsystem, Colon()))
		else
			push!(subsystems, model_or_subsystem)
		end
	end
	return CompositeModel(subsystems...)
end

get_friction_models(system::Vector{<:Subsystem}) = @view system[findall(x->isa(x.model, FrictionModels.ElectronicFrictionProvider), system)]
get_friction_models(system::CompositeModel) = get_friction_models(system.subsystems)
get_pes_models(system::Vector{<:Subsystem}) = @view system[findall(x->isa(x.model, ClassicalModels.ClassicalModel) || isa(x.model, QuantumModels.QuantumModel), system)]
get_pes_models(system::CompositeModel) = get_pes_models(system.subsystems)

dofs(system::CompositeModel) = 1:system.ndofs
ndofs(system::CompositeModel) = system.ndofs

"""
Subsystem combination logic - We only want to allow combination of subsystems:

? 1. with the same number of degrees of freedom

2. without overlapping indices (Build a different type of CompositeModel) to handle these cases separately. 
"""
function check_models(subsystems::Subsystem...)
	systems=vcat(subsystems...)
	# Check unique assignment of potential and derivative to each atom index
	pes_models = get_pes_models(systems)
	pes_model_indices = vcat([subsystem.indices for subsystem in pes_models]...)
	if length(unique(pes_model_indices)) != length(pes_model_indices)
		error("Overlapping indices detected for the assignment of potential energy surfaces.")
	end
	# Check for unique assignment of friction to each atom index
	model_has_friction=systems[findall(x->isa(x.model, FrictionModels.ElectronicFrictionProvider), systems)]
	friction_indices = vcat([subsystem.indices for subsystem in model_has_friction]...)
	if length(unique(friction_indices)) != length(friction_indices)
		error("Overlapping indices detected for the assignment of friction models.")
	end
	
	# Check for different numbers of degrees of freedom in each subsystem
	dofs = [ndofs(subsystem.model) for subsystem in subsystems]
	if length(unique(dofs)) != 1
		error("Subsystems must have the same number of degrees of freedom.")
	end
	return systems, unique(dofs)[1]
end

# Subsystem evaluation of model functions
function potential(system::CompositeModel, R::AbstractMatrix)
	pes_models=get_pes_models(system)
	total_potential_energy=potential(pes_models[1], R)
	for subsystem in pes_models[2:end]
		total_potential_energy+=potential(subsystem, R)
	end
	return total_potential_energy
end

function potential!(system::CompositeModel, V::Real, R::AbstractMatrix)
	return V = potential(system::CompositeModel, R::AbstractMatrix)
end

function derivative!(system::CompositeModel, D::AbstractMatrix, R::AbstractMatrix)
	for subsystem in get_pes_models(system)
		@debug "Accessing D[$(dofs(subsystem)), $(subsystem.indices)]"
		subsystem_derivative = derivative(subsystem, R)
		D[dofs(subsystem), subsystem.indices] .= subsystem_derivative[dofs(subsystem), subsystem.indices]
	end
end

function derivative(system::CompositeModel, R::AbstractMatrix)
	total_derivative=zero(R)
	derivative!(system, total_derivative, R)
	return total_derivative
end

function FrictionModels.friction!(system::CompositeModel, F::AbstractMatrix, R::AbstractMatrix)
	for subsystem in get_friction_models(system)
		if subsystem.indices == Colon()
			eft_indices=vcat([[(j-1)*ndofs(subsystem.model)+i for i in dofs(subsystem.model)] for j in 1:size(R,2)]...) # Size of friction tensor from positions if applying friction to entire system. 
		else
			eft_indices=vcat([[(j-1)*ndofs(subsystem.model)+i for i in dofs(subsystem.model)] for j in subsystem.indices]...)
		end
		FrictionModels.friction!(subsystem, view(F, eft_indices, eft_indices), R)
	end
end

function FrictionModels.friction(system::CompositeModel, R::AbstractMatrix)
	F=zeros(eltype(R), length(R), length(R))
	FrictionModels.friction!(system, F, R)
	return F
end
	
# Overload friction function to map friction_atoms to Subsystem indices
function FrictionModels.friction!(system::Subsystem{<:FrictionModels.DiagonalFriction}, F::AbstractMatrix, R::AbstractMatrix)
	F .= get_friction_matrix(system.model, R)
end

# Overload friction function to map friction_atoms to Subsystem indices
function FrictionModels.friction!(system::Subsystem{<:FrictionModels.TensorialFriction}, F::AbstractMatrix, R::AbstractMatrix)
	indices = FrictionModels.friction_matrix_indices(collect(1:length(system.model.friction_atoms)), NQCModels.ndofs(system))
	F[indices, indices] .= get_friction_matrix(system.model, R)
end

