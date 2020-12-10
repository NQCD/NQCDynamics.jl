export DiffusionAndersonHolstein
export ScatteringAndersonHolstein
export IESH_AnHol_Hamiltonian
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

        #V0 = get_V0
        #D0 = get_D0
        #V = get_potential
        #D = get_derivative
        new(N+1, V0, zero, V, D)
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

        #V0 = get_V0
        #D0 = get_D0
        #V = get_potential
        #D = get_derivative
        new(N+1, V0, zero, V, D)
    end
end


@doc raw"""
    IESH_AnHol_Hamiltonian

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
struct IESH_AnHol_Hamiltonian <: AnalyticModel

    @add_standard_fields

    function IESH_AnHol_Hamiltonian(
        ;N=100::Int64,
        D_0=1.0,
        Bg=1.0,
        Cg=1.0,
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

        #function eslab_states(
        #    c;
        #    E_range = 7.0, 
        #    nsteps = 40
        #    )
        #    # Creates distribution of evenly spaced eigenstates.
        #    e_state = -E_range/2. + (E_range/(nsteps-1))*(c-1)
        #    return e_state
        #
        #end # function eslab_states

        function get_potential(q)::Hermitian

            get_potential = zeros(N+2, N+2)

            get_potential[1,1] = V0(q)
            get_potential[2,2] = V1(q)
            γ_q = γ2(q)
            k = 0
            for i=3:N+1
                get_potential[i,i] = alpha
                # To set the diagonal elements manually
                #k = k + 1
                #V[i,i] = eslab_states(k, E_range=E_range, nsteps=N)
                get_potential[i,i+1] = beta
                get_potential[1,i] = 0.0 # For now, no coupling of GS w/ electrons
                get_potential[2,i] = γ_q
            end
            get_potential[end,end] = alpha
            #V[end,end] = eslab_states(k, E_range=E_range, nsteps=N)
            get_potential[1,end] = 0.0 # For now, no coupling of GS w/ electrons
            get_potential[2,end] = γ_q
            get_potential[end,2] = beta
            get_potential[2,end] = beta

            return Hermitian(get_potential)
        end # End function to calculate energy Hamiltonian V

        #D(q) = Hermitian(zeros(N+1, N+1))
        function get_derivative(q)::Hermitian

            get_derivative = zeros(N+2, N+2)

            get_derivative[1,1] = dV0dr(q)
            get_derivative[2,2] = dV1dr(q)
            γ_q = dγ2dr(q)
            for i=3:N+1
                V[i,i] = 0.0 # dalpha/dr = 0
                V[i,i+1] = 0.0 # dbeta/dr = 0
                V[1,i] = 0.0 # For now, no coupling of GS w/ electrons
                V[2,i] = γ_q
            end
            get_derivative[end,end] = 0.0 # dalpha/dr = 0
            get_derivative[1,end] = 0.0 # For now, no coupling of GS w/ electrons
            get_derivative[2,end] = γ_q
            get_derivative[end,2] = 0.0 # dbeta/dr = 0
            get_derivative[2,end] = 0.0 # dbeta/dr = 0

            return Hermitian(get_derivative)
        end # End function to calculate energy Hamiltonian V

        # This is probably a bad idea, since it assumes the GS PES
        get_D0 = dV0dr
        get_V0 = V0
        # get_V0 is the energy of the GS,
        # get_D0 is the force (position derivative) of the GS
        # get_potential is the diabatic Hamiltonian
        # get_derivative is the diabatic Force matrix
        new(N+2, get_V0, get_D0, get_potential, get_derivative)
    end # function IESH_AnHol_Hamiltonian
end # struc IESH_AnHol_Hamiltonian