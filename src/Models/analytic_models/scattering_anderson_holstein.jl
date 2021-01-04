export ScatteringAndersonHolstein

@doc raw"""
    ScatteringAndersonHolstein

    Anderson-Holstein Hamiltonian for IESH.

```math
V_0(R) = D_e \left(e^{-2aR} - 2e^{-aR}\right)
```
```math
V_1(R) = V_\mathrm{image}(q) + A_1*\left[exp\left(-a_1*q\right)-exp\left(-a_1*q_\mathrm{cut}\right)\right]
```
```math
V_\mathrm(q) = - D_i/\left(C_i^2 + q^2\right)
```
```math
γ(R) = \frac{B}{1 + \exp\left[\frac{r+R}{C}\right]}
```
"""
struct ScatteringAndersonHolstein <: DiabaticModel

    n_states::UInt8
    potential!::Function
    derivative!::Function

    function ScatteringAndersonHolstein(
        ;N::Integer=100,
        D_0=1.0,
        Bg=1.0,
        Cg=0.1,
        shift=1.0,
        a0=1.0,
        D_1 = 1.0,
        C_1 = 1.0,
        A_1 = 1.0,
        a1 = 1.0,
        qcut = 10.0,
        E_range = 7.0,
        alpha=1,
        beta=1
        )
        
        # U0(q) ground state PES, where q is the distance
        V0(q) = D_0*(exp(-2a0*q) - 2*exp(-a0*q))
        # Its derivative:
        dV0dr(q) = -2*D_0*a0*(exp(-a0*q)-exp(-2*a0*q))
        
        # Excited state PES, based on Shenvi, Tully JCP 2009
        # first image charge and then excited state
        # With derivatives
        V_image(q) = - D_1/(C_1^2 + q^2)
        dV_imagedr(q) = -(C_1^2+q^2+2*D_1*q)/(C_1^2 + q^2)^2
        V1(q) = V_image(q) + A_1*(exp(-a1*q)-exp(-a1*qcut))
        dV1dr(q) = -1.0*(V_imagedr(q) -a1*A_1*exp(-a1*q))
        
        # Coupling between particle and slab states
        # With derivatives
        γ2(q)  = Bg/(1+exp((shift+q)/Cg))
        dγ2dr(q) = -(Bg*exp((shift+q)/Cg))/(Cg*(1+exp((shift+q)/Cg))^2)

        """
        Diabatic potential matrix
        """
        function potential!(V::Hermitian, R::AbstractMatrix)
            V[1,1] = V0(R[1])
            V[2,2] = V1(R[1])
            γ_q = γ2(R[1])
            for i=3:N+1
                V[i,i] = alpha
                V.data[i,i+1] = beta
                V.data[2,i] = γ_q
            end
            V[end,end] = alpha
            V.data[2,end] = γ_q
            V.data[3,end] = beta
        end

        """
        Elementwise position derivative of the above `potential!`
        """
        function derivative!(D::AbstractMatrix{<:Hermitian}, R::AbstractArray)
            D[1,1] = dV0dr(R[1])
            D[2,2] = dV1dr(R[1])
            γ_q = dγ2dr(R1[1])
            for i=3:N+2
                D.data[2,i] = γ_q
            end
        end

        new(N+2, potential!, derivative!)
    end
end