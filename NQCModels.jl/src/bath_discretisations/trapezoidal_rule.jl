
"""
    TrapezoidalRule{B,T} <: WideBandBathDiscretisation

Discretise wide band continuum using trapezoidal rule.
Leads to evenly spaced states and constant coupling.
"""
struct TrapezoidalRule{B,T} <: WideBandBathDiscretisation
    bathstates::B   # ϵ
    bathcoupling::T # V(ϵ,x̃) 
end

function TrapezoidalRule(M, bandmin, bandmax)
    ΔE = bandmax - bandmin
    bathstates = range(bandmin, bandmax, length=M)
    bathcoupling = sqrt(ΔE / M)
    return TrapezoidalRule(bathstates, bathcoupling)
end

"""
The **`GapTrapezoidalRule`** struct, a subtype of `WideBandBathDiscretisation`, models a wide-band limit bath with an energy gap.

## Fields

- **`bathstates::Vector{T}`**: The discretized energy levels of the bath.
- **`bathcoupling::T`**: The coupling strength between the system and the bath states.

## Constructor

```julia
GapTrapezoidalRule(M, bandmin, bandmax, gapwidth)
```
This constructor validates inputs and initializes a GapTrapezoidalRule object.

## Arguments

- `M::Int`: The total number of bath states; must be even.
- `bandmin::Real`: The minimum energy of the wide band.
- `bandmax::Real`: The maximum energy of the wide band.
- `gapwidth::Real`: The width of the central energy gap. Must be non-negative and smaller than `bandmax - bandmin`.

## Details

If `gapwidth` is zero, `bathstates` are uniformly distributed across the entire band. Otherwise, the band is split, and `M` states are evenly distributed across the two resulting halves, avoiding the gap.

The `bathcoupling` is calculated from the effective bandwidth (`bandmax - bandmin - gapwidth`) and `M`.

"""
struct GapTrapezoidalRule{T} <: WideBandBathDiscretisation
    bathstates::Vector{T}   # ϵ
    bathcoupling::T # V(ϵ,x̃) 
end

function GapTrapezoidalRule(M, bandmin, bandmax, gapwidth)
    gapwidth < bandmax - bandmin || throw(error("The gap width must be smaller than the band width."))
    gapwidth >= 0 || throw(error("The gap width must be positive."))
    M % 2 == 0 || throw(error("The number of states `M` must be even."))
    ΔE = bandmax - bandmin
    centre = (bandmax + bandmin) / 2

    if gapwidth == 0
        bathstates = collect(range(bandmin, bandmax, length=M))
    else
        bathstates_left = range(bandmin, centre - gapwidth/2, length=Int(M/2))
        bathstates_right = range(centre + gapwidth/2, bandmax, length=Int(M/2))
        bathstates = vcat(bathstates_left, bathstates_right)
    end

    bathcoupling = sqrt((ΔE-gapwidth) / M)
    return GapTrapezoidalRule(bathstates, bathcoupling)
end
