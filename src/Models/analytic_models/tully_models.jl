export TullyModelOne
export TullyModelTwo

"""
The 2-state Tully model 1
"""
struct TullyModelOne <: DiabaticModel

    n_states::UInt8
    potential!::Function
    derivative!::Function

    function TullyModelOne(A=0.01, B=1.6, C=0.005, D=1.0)

        function potential!(V::Hermitian, R::AbstractMatrix)
            q = R[1]
            if q > 0
                V[1,1] = A * (1 - exp(-B*q))
            else
                V[1,1] = -A * (1 - exp(B*q))
            end
            V[2,2] = -V[1,1]
            V.data[1,2] = C * exp(-D*q^2) # sets both off-diagonals as V is Hermitian
        end

        function derivative!(derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
            q = R[1]
            derivative[1][1,1] = A * B * exp(-B * abs(q))
            derivative[1][2,2] = -D[1][1,1]
            derivative[1].data[1,2] = -2 * C * D * q * exp(-D*q^2)
        end

        new(2, potential!, derivative!)
    end
end

struct TullyModelTwo <: DiabaticModel

    n_states::UInt8
    potential!::Function
    derivative!::Function

    function TullyModelTwo(A=0.1, B=0.28, C=0.015, D=0.06, E=0.05)

        function potential!(V::Hermitian, R::AbstractMatrix)
            q = R[1]
            V[1,1] = 0
            V[2,2] = -A*exp(-B*q^2) + E
            V.data[1,2]= C * exp(-D*q^2)
        end

        function derivative!(derivative::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
            q = R[1]
            derivative[1][1,1] = 0
            derivative[1][2,2] = 2*A*B*q*exp(-B*q^2)
            derivative[1].data[1,2] = -2*C*D*q*exp(-D*q^2)
        end

        new(2, potential!, derivative!)
    end
end
