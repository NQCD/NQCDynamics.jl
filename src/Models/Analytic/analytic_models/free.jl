export Free

"""
A model that has zero external potential everywhere.
"""
struct Free <: AnalyticModel

    @add_standard_fields

    Free() = new(1, x->0, x->0, zero_hermitian, zero_hermitian)
end
