"""
    ShenviGaussLegendre{T}

Defined as described by Shenvi et al. in J. Chem. Phys. 130, 174107 (2009).
Use Gauss-Legendre quadrature to discretise the bath continuum in two halves across the given bandwidth, split at the Fermi level.
The position of the negative sign for the state energy level has been moved to ensure the states are sorted from lowest to highest.
"""
struct ShenviGaussLegendre{T} <: BathDiscretisationScheme
    bathstates::Vector{T}
    bathcoupling::Vector{T}
    discretisationtype::Symbol
end

"""
    ShenviGaussLegendre(M::Int64, bandmin, bandmax)

Args:
- `M`: number of states in discretisation  
- `bandmin`: minimum of discretised energy range
- `bandmax`: maximum of discretised energy range
"""
function ShenviGaussLegendre(M::Int64, bandmin, bandmax)
    M % 2 == 0 || throw(error("The number of states `M` must be even."))
    knots, weights = gausslegendre(div(M, 2))
    centre = (bandmax + bandmin) / 2

    bathstates = zeros(M)
    for i in eachindex(knots)
        bathstates[i] = (centre - bandmin)/2 * knots[i] + (bandmin + centre) / 2
    end
    for i in eachindex(knots)
        bathstates[i+length(knots)] = (bandmax - centre)/2 * knots[i] + (bandmax + centre) / 2
    end

    bathcoupling = zeros(M)
    for (i, w) in enumerate(weights)
        bathcoupling[i] = sqrt((centre - bandmin)/2  * w)
    end
    for (i, w) in enumerate(weights)
        bathcoupling[i+length(weights)] = sqrt((bandmax - centre)/2  * w)
    end

    discretisationtype = :ShenviGaussLegendre

    return ShenviGaussLegendre(bathstates, bathcoupling, discretisationtype)
end

"""
    ReferenceGaussLegendre{T}

Implementation translated from Fortran code used for simulations of Shenvi et al. in J. Chem. Phys. 130, 174107 (2009).
Two differences from ShenviGaussLegendre:
- Position of minus sign in energy levels has been corrected.
- Division by sqrt(ΔE) in the coupling. 
"""
struct ReferenceGaussLegendre{T} <: BathDiscretisationScheme
    bathstates::Vector{T}
    bathcoupling::Vector{T}
end

function ReferenceGaussLegendre(M, bandmin, bandmax)
    M % 2 == 0 || throw(error("The number of states `M` must be even."))
    @warn "(ReferenceGaussLegendre) The division by sqrt(ΔE) in the coupling is incorrect. Use ShenviGaussLegendre instead."
    knots, weights = gausslegendre(div(M, 2))
    ΔE = bandmax - bandmin

    bathstates = zeros(M)
    for i in eachindex(knots)
        bathstates[i] = ΔE/2 * (-1/2 + knots[i]/2)
    end
    for i in eachindex(knots)
        bathstates[i+length(knots)] = ΔE/2 * (1/2 + knots[i]/2)
    end

    bathcoupling = zeros(M)
    for (i, w) in zip(eachindex(bathcoupling), Iterators.cycle(weights))
        bathcoupling[i] = sqrt(w) / 2
    end

    return ReferenceGaussLegendre(bathstates, bathcoupling)
end

