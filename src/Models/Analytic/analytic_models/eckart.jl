export Eckart

"""
Symmetric Eckart barrier.
"""
struct Eckart <: AnalyticModel

    @add_standard_fields

    function Eckart(v0, a)

        V0(q)::Float64 = v0 / cosh(q/a)^2
        D0(q)::Float64 = -2 * v0 / a * tanh(q/a) * sech(q/a)^2

        new(V0, D0, zero_hermitian, zero_hermitian)
    end
end
