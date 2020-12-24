export DoubleWell

"""
The 2-state double well model often used in the literature.
"""
struct DoubleWell <: DiabaticModel

    n_states::UInt8
    potential!::Function
    derivative!::Function

    function DoubleWell(;M=1, ω=1, γ=1, Δ=1)

        V0(R) = 0.5 * M * ω^2 * R^2
        D0(R) = M * ω^2 * R
        function potential!(V::Hermitian, R::AbstractMatrix)
            v = sqrt(2)*γ*R[1] # Access only R[1] as this is a 1D model
            V[1,1] = v
            V[2,2] = -v
            V.data[1,2] = Δ/2 # Sets both off-diagonals as V is Hermitian
            V += I*V0(R[1])
        end

        function derivative!(D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
            v = sqrt(2)*γ
            D[1][1,1] = v
            D[1][2,2] = -v
            D[1] += I*D0(R[1])
        end

        new(2, potential!, derivative!)
    end
end
