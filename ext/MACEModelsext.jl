module MACEModelsext
"""
To take advantage of batch evaluation in ring-polymer MD simulations, we need to overload potential
and derivative evaluation in the ring polymer calculators in NQCDynamics. 
"""

using MACEModels
using NQCDynamics
using NQCCalculators

function NQCCalculators.update_potential!(cache::RingPolymer_ClassicalModel_Cache{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    potential_vector = @views [R[:, :, i] for i in axes(R, 3)]
    cache.potential = NQCModels.potential(cache.model, cache.model.atoms, potential_vector, cache.model.cell)
    return nothing
end

function NQCCalculators.update_derivative!(cache::RingPolymer_ClassicalModel_Cache{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    derivative_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.derivative!(cache.model, cache.model.atoms, cache.derivative, derivative_vector, cache.model.cell)
    return nothing
end

function NQCCalculators.update_potential!(cache::RingPolymer_ClassicalFrictionModel_Cache{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    @debug "RPMDEF accelerated potential evaluation"
    potential_vector = @views [R[:, :, i] for i in axes(R, 3)]
    cache.potential = NQCModels.potential(cache.model, cache.model.atoms, potential_vector, cache.model.cell)
    return nothing
end

function NQCCalculators.update_derivative!(cache::RingPolymer_ClassicalFrictionModel_Cache{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    @debug "RPMDEF accelerated derivative evaluation"
    derivative_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.derivative!(cache.model, cache.model.atoms, cache.derivative, derivative_vector, cache.model.cell)
    return nothing
end

end
