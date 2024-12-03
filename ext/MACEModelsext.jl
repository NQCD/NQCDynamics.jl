module MACEModelsext
"""
To take advantage of batch evaluation in ring-polymer MD simulations, we need to overload potential
and derivative evaluation in the ring polymer calculators in NQCDynamics. 
"""

using MACEModels
using NQCDynamics: Calculators

function Calculators.evaluate_potential!(calc::RingPolymerAdiabaticCalculator{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    calc.stats[:potential] += 1
    potential_vector = @views [R[:, :, i] for i in axes(R, 3)]
    calc.potential = NQCModels.potential(calc.model, calc.model.atoms, potential_vector, calc.model.cell)
    return nothing
end

function Calculators.evaluate_derivative!(calc::RingPolymerAdiabaticCalculator{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    calc.stats[:derivative] += 1
    derivative_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.derivative!(calc.model, calc.model.atoms, calc.derivative, derivative_vector, calc.model.cell)
    return nothing
end

function Calculators.evaluate_potential!(calc::RingPolymerFrictionCalculator{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    calc.stats[:potential] += 1
    potential_vector = @views [R[:, :, i] for i in axes(R, 3)]
    calc.potential = NQCModels.potential(calc.model, calc.model.atoms, potential_vector, calc.model.cell)
    return nothing
end

function Calculators.evaluate_derivative!(calc::RingPolymerFrictionCalculator{T,MACEModel}, R::AbstractArray{S,3}) where {T,S}
    calc.stats[:derivative] += 1
    derivative_vector = @views [R[:, :, i] for i in axes(R, 3)]
    NQCModels.derivative!(calc.model, calc.model.atoms, calc.derivative, derivative_vector, calc.model.cell)
    return nothing
end

	
end