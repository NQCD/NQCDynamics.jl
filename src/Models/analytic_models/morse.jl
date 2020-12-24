export MorsePotential

"""
The 3-state morse potential J.Chem.Phys. 150, 244102(2019)
"""
struct MorsePotential <: AnalyticModel

    @add_standard_fields

    function MorsePotential()

        d1 = 0.02; d2 = 0.02; d3 = 0.003
        α1 = 0.4; α2 = 0.65; α3 = 0.65
        r1 = 4.0; r2 = 4.5; r3 = 6.0
        c1 = 0.02; c2 = 0.0; c3 = 0.02

        a12 = 0.005; a13 = 0.005; a23 = 0.0
        α12 = 32.0; α13 = 32.0; α23 = 0.0
        r12 = 3.40; r13 = 4.97; r23 = 0.0

        function V_ii(x, d, α, r, c)
            return d * (1 - exp(-α*(x-r)))^2 + c
        end

        function V_ij(x, a, α, r)
            return a * exp(-α*(x-r)^2)
        end

        function morse_potential(q)::Hermitian

            v1 = V_ii(q, d1, α1, r1, c1)
            v2 = V_ii(q, d2, α2, r2, c2)
            v3 = V_ii(q, d3, α3, r3, c3)

            v12 = V_ij(q, a12, α12, r12)
            v13 = V_ij(q, a13, α13, r13)
            v23 = V_ij(q, a23, α23, r23)

            return Hermitian([v1 v12 v13; v12 v2 v23; v13 v23 v3])
        end

        function D_ii(x, d, α, r)
            ex = exp(-α*(x-r))
            return 2 * d * α * (ex - ex^2)
        end

        function D_ij(x, a, α, r)
            return -2 * a * α * (x-r) * exp(-α*(x-r)^2)
        end

        function morse_derivative(q)::Hermitian

            v1 = D_ii(q, d1, α1, r1)
            v2 = D_ii(q, d2, α2, r2)
            v3 = D_ii(q, d3, α3, r3)

            v12 = D_ij(q, a12, α12, r12)
            v13 = D_ij(q, a13, α13, r13)
            v23 = D_ij(q, a23, α23, r23)

            return Hermitian([v1 v12 v13; v12 v2 v23; v13 v23 v3])
        end

        new(3, zero, zero, morse_potential, morse_derivative)
    end
end
