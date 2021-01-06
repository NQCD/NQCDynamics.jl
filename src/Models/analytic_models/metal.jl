export DiffusionMetal
export ScatteringMetal
using LinearAlgebra

@doc raw"""
    DiffusionMetal <: Model

A simple diffusive variety of the Anderson-Holstein Hamiltonian.

```math
U_0 = A\cos\left(\frac{2π}{a} (x-a/2)\right)\\
V_i(R) = B\cos^2\left(\frac{π}{a}x\right)\\
\epsilon_d(R) = 0
```
The electronic levels are evenly spaced over a given region.   
"""
struct DiffusionMetal <: AnalyticModel

    @add_standard_fields

    function DiffusionMetal(
        ;N=100::Integer,
        band_range::Tuple=(-2, 2),
        a=1,
        B=1)

        U0(x) = cos(2π/a*(x-a/2))
        D0(x) = -sin(2π/a*(x-a/2)) * 2π/a

        h(x) = 0
        dh(x) = 0

        coupling(x) = B*cos(π/a*x)^2
        dcoupling(x) = B*2*cos(π/a*x)*-sin(π/a*x)*π/a

        energies = collect(range(band_range..., length=N))

        function V(q)::Hermitian

            v = zeros(Float64, N+1, N+1)
            v[1,1] = h(q)
            v[2:end,2:end] = diagm(energies)
            v[1,2:end] .= coupling(q) / sqrt(N)

            return Hermitian(v)
        end

        function D(q)::Hermitian

            v = zeros(Float64, N+1, N+1)
            v[1,1] = dh(q)
            v[1,2:end] .= dcoupling(q) / sqrt(N)

            return Hermitian(v)
        end

        new(N+1, U0, D0, V, D)
    end
end

@doc raw"""
    ScatteringMetal <: Model

A simple scattering variety of the Anderson-Holstein Hamiltonian.

```math
U_0(R) = D_e \left(\exp(-2aR) -2\exp(-aR) \right)\\
\epsilon_d(R) = h \exp(-ax) - U_0(R)\\
V_i(R) = \frac{B}{1+\exp(\frac{δ + R}{C})}
```
The electronic levels are evenly spaced over a given region.   
"""
struct ScatteringMetal <: AnalyticModel

    @add_standard_fields

    function ScatteringMetal(
        ;N=100::Int,
        band_range::Tuple=(-2,2),
        a=1,
        B=0.5,
        C=0.5,
        D_e=1,
        height=1,
        shift=-1)

        U0(q) = D_e*(exp(-2a*q) - 2exp(-a*q))
        D0(q) = D_e*(-2a*exp(-2a*q) + 2a*exp(-a*q))

        h(x) = height*exp(-a*x) - U0(x)
        dh(x) = -a*height*exp(-a*x) - D0(x)

        coupling(q)  = B/(1+exp((shift+q)/C))
        # Correct this function please sir 
        dcoupling(q)  = B/(1+exp((shift+q)/C))

        energies = collect(range(band_range..., length=N))

        function V(q)::Hermitian

            v = zeros(Float64, N+1, N+1)
            v[1,1] = h(q)
            v[2:end,2:end] = diagm(energies)
            v[1,2:end] .= coupling(q) / sqrt(N)

            return Hermitian(v)
        end

        function D(q)::Hermitian

            d = zeros(Float64, N+1, N+1)
            d[1,1] = dh(q)
            d[1,2:end] .= dcoupling(q) / sqrt(N)

            return Hermitian(d)
        end

        new(N+1, U0, D0, V, D)
    end
end
