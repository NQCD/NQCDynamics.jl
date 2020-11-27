export Free

"""
A model that has zero external potential everywhere.
"""
struct Free <: AnalyticModel

    @add_standard_fields

    function Free()
        new(x->0, x->0, zero_hermitian, zero_hermitian)
    end
end
