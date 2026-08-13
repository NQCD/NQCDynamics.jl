"""
    FullGaussLegendre{T} <: WideBandBathDiscretisation

Use Gauss-Legendre quadrature to discretise the continuum across the entire band width.
This is similar to the ShenviGaussLegendre except that splits the continuum at the Fermi level into two halves.
"""
struct FullGaussLegendre{T} <: WideBandBathDiscretisation
    bathstates::Vector{T}
    bathcoupling::Vector{T}
end

function FullGaussLegendre(M, bandmin, bandmax)

    knots, weights = gausslegendre(M)
    ΔE = bandmax - bandmin

    bathstates = zeros(M)
    for i in eachindex(knots, bathstates)
        bathstates[i] = ΔE/2 * knots[i]
    end

    bathcoupling = zeros(M)
    for i in eachindex(bathcoupling, weights)
        bathcoupling[i] = sqrt(weights[i] * ΔE/2)
    end

    return FullGaussLegendre(bathstates, bathcoupling)
end



"""
The **`GapGaussLegendre`** struct, a subtype of `WideBandBathDiscretisation`, models a wide-band limit bath with a central energy gap using Gauss–Legendre quadrature for discretization.

## Fields

- **`bathstates::Vector{T}`**: The discretized energy levels of the bath.
- **`bathcoupling::Vector{T}`**: The coupling strengths to the bath states, derived from Gauss–Legendre weights.

## Constructor

```julia
GapGaussLegendre(M, bandmin, bandmax, gapwidth)
```

This constructor validates input arguments and initializes a `GapGaussLegendre` object using Gauss–Legendre quadrature for state placement and coupling computation.

## Arguments

- `M::Int`: The total number of bath states; must be even.
- `bandmin::Real`: The minimum energy of the wide band.
- `bandmax::Real`: The maximum energy of the wide band.
- `gapwidth::Real`: The width of the central energy gap. Must be non-negative and smaller than `bandmax - bandmin`.

## Details

The energy band is split into two symmetric regions around the gap, and `M/2` states are placed in each region using Gauss–Legendre nodes mapped to the corresponding subinterval. The gap itself is avoided.

The `bathcoupling` vector is computed from the square roots of the Gauss–Legendre weights, scaled by the effective half-bandwidth. This ensures accurate integration over the energy band while excluding the gap.
"""
struct GapGaussLegendre{T} <: WideBandBathDiscretisation
    bathstates::Vector{T}
    bathcoupling::Vector{T}
end

function GapGaussLegendre(M, bandmin, bandmax, gapwidth)
    gapwidth < bandmax - bandmin || throw(error("The gap width must be smaller than the band width."))
    gapwidth >= 0 || throw(error("The gap width must be positive."))
    M % 2 == 0 || throw(error("The number of states `M` must be even."))
    knots, weights = gausslegendre(div(M, 2))
    centre = (bandmax + bandmin) / 2
    ΔE = bandmax - bandmin
    bandlengths = (ΔE - gapwidth)/2

    bathstates = zeros(M)
    for i in eachindex(knots)
        bathstates[length(knots)+1-i] = ((ΔE-gapwidth)/2 * (1 + knots[i]) + gapwidth) * (-1/2) + centre
    end
    for i in eachindex(knots)
        bathstates[i+length(knots)] = ((ΔE-gapwidth)/2 * (1 + knots[i]) + gapwidth) * 1/2 + centre
    end

    # coupling is the sqrt() of the weights times the bandlengths [valence band's coupling ; conduction band's coupling]
    bathcoupling = [sqrt.(weights.*bandlengths/2); sqrt.(weights.*bandlengths/2)]

    return GapGaussLegendre(bathstates, bathcoupling)
end
