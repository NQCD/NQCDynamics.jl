export Harmonic

"""
The 1-state harmonic model
"""
struct Harmonic <: AnalyticModel

    @add_standard_fields

    function Harmonic(mass, omega, r_0)

        V(q) = 0.5 * mass * omega ^ 2 * (q - r_0) ^2
        D(q) = mass * omega ^ 2 * (q - r_0)

        new(1, V, D, zero_hermitian, zero_hermitian)
    end
end
