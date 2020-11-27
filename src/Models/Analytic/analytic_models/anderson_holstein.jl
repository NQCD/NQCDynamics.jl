export DiffusionAndersonHolstein
export ScatteringAndersonHolstein
using LinearAlgebra

@doc raw"""
    DiffusionAndersonHolstein

Diffusive Anderson-Holstein Hamiltonian with a tight-binding electronic bath.

```math
U_0(R) = A\cos[\frac{2π}{a} (R-\frac{a}{2})]
```
```math
γ_i(R) = B\exp\left[-\frac{(r_i - R)^2}{C}\right]
```
```math
\epsilon_d(R) = 0
```

The values of $r_i$ are equally spaced by $a$, leading to a local coupling to
the local basis functions on each metallic atom.
"""
struct DiffusionAndersonHolstein <: AnalyticModel

    @add_standard_fields

    function DiffusionAndersonHolstein(
        ;N=100::Int64,
        α=0,
        β=1,
        a=1,
        A=1,
        B=1,
        C=1)

        V0(q) = A*cos(2*π/a * (q-a/2))

        ϵd(q) = 0
        γ(q, n) = B*exp(-(n*a-q)^2/C)

        function V(q)::Hermitian

            V = zeros(N+1, N+1)

            V[1,1] = ϵd(q)
            for i=2:N
                V[i,i] = α
                V[i,i+1] = β
                V[1,i] = γ(q, i-2-N÷2)
            end
            V[end,end] = α
            V[1,end] = γ(q, N-1-N÷2)
            V[end,2] = β
            V[2,end] = β

            return Hermitian(V)
        end

        D(q) = Hermitian(zeros(N+1, N+1))

        new(V0, zero, V, D)
    end
end


@doc raw"""
    ScatteringAndersonHolstein

Scattering Anderson-Holstein Hamiltonian with a tight-binding electronic bath.

```math
U_0(R) = D_e \left(e^{-2aR} - 2e^{-aR}\right)
```
```math
\epsilon_d(R) = h\exp(-ax) - V_0(R)
```
```math
γ(R) = \frac{B}{1 + \exp\left[\frac{r+R}{C}\right]}
```
"""
struct ScatteringAndersonHolstein <: AnalyticModel

    @add_standard_fields

    function ScatteringAndersonHolstein(
        ;N=100::Int64,
        α=1,
        β=1,
        a=1,
        D_e=1,
        B=1,
        C=1,
        shift=1,
        height=1
        )

        V0(q) = D_e*(exp(-2a*q) - 2exp(-a*q))

        ϵd(q) = height*exp(-a*q) - V0(q)
        γ(q)  = B/(1+exp((shift+q)/C))

        function V(q)::Hermitian

            V = zeros(N+1, N+1)

            V[1,1] = ϵd(q)
            γ_q = γ(q)
            for i=2:N
                V[i,i] = α
                V[i,i+1] = β
                V[1,i] = γ_q
            end
            V[end,end] = α
            V[1,end] = γ_q
            V[end,2] = β
            V[2,end] = β

            return Hermitian(V)
        end

        D(q) = Hermitian(zeros(N+1, N+1))

        new(V0, zero, V, D)
    end
end