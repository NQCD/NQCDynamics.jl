module Analytic

using LinearAlgebra: Hermitian
using ..Models

abstract type AnalyticModel <: Models.Model end

"""
    @add_standard_fields

Add the fields required for all of the models.

This should be used for every model to ensure they have all the necessary fields.
The four fields are functions that return the state-independent potential
and the state-dependent diabatic matrix along with their associated derivatives.
"""
macro add_standard_fields()
    esc(quote
        get_V0::Function
        get_D0::Function
        get_potential::Function
        get_derivative::Function
    end)
end

"""
    zero_hermitian(A)

Return a Hermitian 1x1 matrix.

If the model is completely state-independent, this function can be used
for the state-dependent potential and derivative.
"""
zero_hermitian(A) = Hermitian(zeros(1,1))

include("analytic_models/double_well.jl")
include("analytic_models/free.jl")
include("analytic_models/harmonic.jl")
include("analytic_models/morse.jl")
include("analytic_models/tully_models.jl")
include("analytic_models/eckart.jl")
include("analytic_models/metal.jl")
include("analytic_models/metallic_chain.jl")
include("analytic_models/anderson_holstein.jl")
include("analytic_models/spin_boson.jl")

include("plot.jl")

end # module