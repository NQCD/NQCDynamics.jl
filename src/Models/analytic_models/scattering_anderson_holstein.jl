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
    N::UInt8
    D_0::Float64
    Bg::Float64
    Cg::Float64
    shift::Float64
    a0::Float64
    D_1::Float64
    C_1::Float64
    A_1::Float64
    a1::Float64
    qcut::Float64
    E_range::Float64
    alpha::Float64
    beta::Float64
end

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
    ScatteringAndersonHolstein(N+2, N, D_0, Bg, Cg, shift, a0, D_1 , C_1, A_1, a1,
        qcut, E_range, alpha, beta)
end

"""
Diabatic potential matrix
"""
function potential!(model::ScatteringAndersonHolstein, V::Hermitian, R::AbstractMatrix)

    # U0(q) ground state PES, where q is the distance
    V0(q) = model.D_0*(exp(-2model.a0*q) - 2*exp(-model.a0*q))
    # Excited state PES, based on Shenvi, Tully JCP 2009
    # first image charge and then excited state
    V_image(q) = -model.D_1/(model.C_1^2 + q^2)
    V1(q) = V_image(q) + model.A_1*(exp(-model.a1*q)-exp(-model.a1*model.qcut))
    # Coupling between particle and slab states
    γ2(q) = model.Bg/(1+exp((model.shift+q)/model.Cg))

    V[1,1] = V0(R[1])
    V[2,2] = V1(R[1])
    γ_q = γ2(R[1])
    for i=3:model.N+1
        V[i,i] = model.alpha
        V.data[i,i+1] = model.beta
        V.data[2,i] = γ_q
    end
    V[end,end] = model.alpha
    V.data[2,end] = γ_q
    V.data[3,end] = model.beta
end

"""
Elementwise position derivative of the above `potential!`
"""
function derivative!(model::ScatteringAndersonHolstein, D::AbstractMatrix{<:Hermitian}, R::AbstractArray)
    dV0dr(q) = -2*model.D_0*model.a0*(exp(-model.a0*q)-exp(-2*model.a0*q))
    dV_imagedr(q) = -(model.C_1^2+q^2+2*model.D_1*q)/(model.C_1^2 + q^2)^2
    dV1dr(q) = -1.0*(dV_imagedr(q) -model.a1*model.A_1*exp(-model.a1*q))
    dγ2dr(q) = -(model.Bg*exp((model.shift+q)/model.Cg))/(model.Cg*(1+exp((model.shift+q)/model.Cg))^2)

    D[1][1,1] = dV0dr(R[1])
    D[1][2,2] = dV1dr(R[1])
    γ_q = dγ2dr(R[1])
    for i=3:model.N+2
        D[1].data[2,i] = γ_q
    end
end
