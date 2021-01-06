export OhmicSpinBoson
export OhmicBosonicBath
export DebyeSpinBoson
export DebyeBosonicBath

"""
Model for bath of bosonic oscillators with Ohmic spectral density coupled to a spin.

Bath is discretised following: 
J. Phys. Chem. A 2003, 107, 2126-2136
and
J. Chem. Phys. 115, 2979 (2001).
"""
struct OhmicSpinBoson <: AnalyticModel

    @add_standard_fields

    """
        OhmicSpinBoson(Nᵇ; ϵ=0, Δ=1, η=0.09, ωᶜ=2.5)
    """
    function OhmicSpinBoson(Nᵇ; ϵ=0, Δ=1, η=0.09, ωᶜ=2.5)

        ωᵐ = 10ωᶜ
        J(ω) = π/2 * η * ω * exp(-ω/ωᶜ)
        ρ(ω) = Nᵇ / ωᶜ * exp(-ω/ωᶜ)/(1 - exp(-ωᵐ/ωᶜ))
        c(ω) = sqrt(2/π * ω * J(ω)/ρ(ω))

        ωⱼ = range(0, ωᵐ, length=Nᵇ+1)[2:end]
        cⱼ = c.(ωⱼ)

        σ_z = Hermitian([1 0; 0 -1])
        σ_x = Hermitian([0 1; 1 0])

        V0(R::Vector)::Number = sum(ωⱼ.^2 .* R.^2)/2
        D0(R::Vector)::Vector = ωⱼ.^2 .* R

        V(R::Vector)::Hermitian = (sum(cⱼ.*R) + ϵ)*σ_z + Δ*σ_x
        D(R::Vector)::Vector{Hermitian} = [j*σ_z for j in cⱼ]

        new(2, V0, D0, V, D)
    end
end

@doc raw"""
The isolated bath that appears in the spin boson model.

```math
H_b = \frac{1}{2}\sum_j p_j^2 + \omega_j^2 q_j^2
```
"""
struct OhmicBosonicBath <: AnalyticModel
    @add_standard_fields

    """
        BosonBath(Nᵇ; ωᶜ=2.5)
    """
    function OhmicBosonicBath(Nᵇ; ωᶜ=2.5)

        ωᵐ = 10ωᶜ
        ωⱼ = range(0, ωᵐ, length=Nᵇ+1)[2:end]

        V0(R::Vector)::Number = sum(ωⱼ.^2 .* R.^2)/2
        D0(R::Vector)::Vector = ωⱼ.^2 .* R
        function D(R::Vector)::Vector
            d = zero_hermitian(R)
            [d for _=1:Nᵇ]
        end

        new(1, V0, D0, zero_hermitian, D)
    end
end

"""
    DebyeSpinBoson <: Model

P. Huo and D. F. Coker, Mol. Phys. 110, 1035 (2012).
Rekik et al. J. Chem. Phys. 138, 144106 (2013).
"""
struct DebyeSpinBoson <: AnalyticModel

    @add_standard_fields

    """
        DebyeSpinBoson(Nᵇ; ϵ=0, Δ=1, η=0.09, ωᶜ=2.5)

    `J(ω) = η * ω * ωᶜ / (ω^2 + ωᶜ^2)`
    """
    function DebyeSpinBoson(Nᵇ; ϵ=0, Δ=1, η=0.09, ωᶜ=2.5)

        ωᵐ = 10ωᶜ
        c(ω) = sqrt(2η * atan(ωᵐ / ωᶜ) / (π * Nᵇ)) * ω
        ω(j) = tan(j * atan(ωᵐ / ωᶜ) / Nᵇ) * ωᶜ

        ωⱼ = ω.(1:Nᵇ)
        cⱼ = c.(ωⱼ)

        σ_z = Hermitian([1 0; 0 -1])
        σ_x = Hermitian([0 1; 1 0])

        V0(R::Vector)::Number = sum(ωⱼ.^2 .* R.^2)/2
        D0(R::Vector)::Vector = ωⱼ.^2 .* R

        V(R::Vector)::Hermitian = (sum(cⱼ.*R) + ϵ)*σ_z + Δ*σ_x
        D(R::Vector)::Vector{Hermitian} = [j*σ_z for j in cⱼ]

        new(2, V0, D0, V, D)
    end
end

"""
    DebyeBosonicBath <: Model

Adiabatic bath to go with the associated spin boson model.

Used for sampling the bath uncoupled to the spin.
"""
struct DebyeBosonicBath <: AnalyticModel

    @add_standard_fields
    ωⱼ::Vector{Real}

    """
        DebyeBosonicBath(Nᵇ; η=0.09, ωᶜ=2.5)
    """
    function DebyeBosonicBath(Nᵇ; η=0.09, ωᶜ=2.5)

        ωᵐ = 10ωᶜ
        ω(j) = tan(j * atan(ωᵐ / ωᶜ) / Nᵇ) * ωᶜ

        ωⱼ = ω.(1:Nᵇ)

        V0(R::Vector)::Number = sum(ωⱼ.^2 .* R.^2)/2
        D0(R::Vector)::Vector = ωⱼ.^2 .* R
        function D(R::Vector)::Vector
            d = zero_hermitian(R)
            [d for _=1:Nᵇ]
        end

        new(1, V0, D0, zero_hermitian, D, ωⱼ)
    end
end