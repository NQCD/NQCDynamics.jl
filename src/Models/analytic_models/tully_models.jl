export TullyModelOne
export TullyModelTwo

"""
The 2-state Tully model 1
"""
struct TullyModelOne <: AnalyticModel

    @add_standard_fields

    function TullyModelOne(A=0.01, B=1.6, C=0.005, D=1.0)

        function V(q)::Hermitian
            if q > 0
                v1 = A * (1 - exp(-B*q))
            else
                v1 = -A * (1 - exp(B*q))
            end
            v2 = -v1
            v12 = C * exp(-D*q^2)
            return Hermitian([v1 v12; v12 v2])
        end

        function Deriv(q)::Hermitian
            v1 = A * B * exp(-B * abs(q))
            v2 = -v1
            v12 = -2 * C * D * q * exp(-D*q^2)
            return Hermitian([v1 v12; v12 v2])
        end

        new(2, zero, zero, V, Deriv)
    end
end

struct TullyModelTwo <: AnalyticModel

    @add_standard_fields

    function TullyModelTwo(A=0.1, B=0.28, C=0.015, D=0.06, E=0.05)

        function V(q)::Hermitian
            v1 = 0
            v2 = -A*exp(-B*q^2) + E
            v12 = C * exp(-D*q^2)
            return Hermitian([v1 v12; v12 v2])
        end

        function Deriv(q)::Hermitian
            v1 = 0
            v2 = 2*A*B*q*exp(-B*q^2)
            v12 = -2*C*D*q*exp(-D*q^2)
            return Hermitian([v1 v12; v12 v2])
        end

        new(2, zero, zero, V, Deriv)
    end
end
