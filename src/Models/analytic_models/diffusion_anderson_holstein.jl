export DiffusionAndersonHolstein

struct DiffusionAndersonHolstein <: DiabaticModel

    n_states::UInt8
    potential!::Function
    derivative!::Function

    function DiffusionAndersonHolstein(
        ;N=100::Integer,
        α=0,
        β=1,
        a=1,
        A=1,
        B=1,
        C=1)

        V0(R) = A*cos(2π/a * R)
        D0(R) = -A*sin(2π/a*R)*2π/a
        γ(R, n) = B*exp(-(n*a-R)^2/C)
        dγ(R, n) = γ(R, n) * 2 * (n*a-R)/C 

        function potential!(V::Hermitian, R::AbstractMatrix)
            V[1,1] = V0(R[1])
            for i=2:N
                V.data[i,i] = α
                V.data[i,i+1] = β
                V.data[1,i] = γ(R[1],i-1-N÷2)
            end
            V.data[1,end] = γ(R[1], N-N÷2)
            V.data[end,end] = α
        end

        function derivative!(D::AbstractMatrix{<:Hermitian}, R::AbstractArray)
            D[1,1] = D0(R[1])
            for i=2:N+1
                D.data[1,i] = dγ(R, n)
            end
        end

        new(N+1, potential!, derivative!)
    end
end