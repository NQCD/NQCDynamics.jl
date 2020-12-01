export DoubleWell

"""
The 2-state double well model often used in the literature.
"""
struct DoubleWell <: AnalyticModel

    @add_standard_fields

    function DoubleWell(mass, omega, delta, k)

        function V(q)::Hermitian
            v0 = 0.5 * mass * omega ^ 2 * q ^2
            v1 = k * q
            v2 = -v1
            return Hermitian([v0+v1 delta; delta v0+v2])
        end

        function D(q)::Hermitian
            d0 = mass * omega ^ 2 * q
            return Hermitian([d0+k 0; 0 d0-k])
        end

        new(2, zero, zero, V, D)
    end

    # function DoubleWell(;m, ω, g, ϵ)
    #
    #     V0(q)::Float64 = 0
    #     D0(q)::Float64 = 0
    #     function double_well_potential(q)::Hermitian
    #         v0 = 0.5*m*ω^2 * q^2
    #         h = ϵ + g*q*sqrt(2*m*ω)
    #         return Hermitian([v0 0; 0 v0+h])
    #     end
    #
    #     function double_well_derivative(q)::Hermitian
    #         d0 = m*ω^2 * q
    #         dh = g*sqrt(2*m*ω)
    #         return Hermitian([d0 0; 0 d0+dh])
    #     end
    #
    #     new(V0, D0, double_well_potential, double_well_derivative,
    #         nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    # end

    """
        DoubleWell(M, ω, γ, Δ)

    """
    function DoubleWell(;M=1, ω=1, γ=1, Δ=1)

        V0(q)::Float64 = 0.5 * M * ω^2 * q^2
        D0(q)::Float64 = M * ω^2 * q
        function V(q)::Hermitian
            v = sqrt(2)*γ*q
            return Hermitian([v Δ/2; Δ/2 -v])
        end

        function D(q)::Hermitian
            v = sqrt(2)*γ
            return Hermitian([v 0; 0 -v])
        end

        new(2, V0, D0, V, D)
    end
end
