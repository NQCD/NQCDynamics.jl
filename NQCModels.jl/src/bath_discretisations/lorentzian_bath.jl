using Unitful, UnitfulAtomic
"""
    lorentzian(energy::Float64, width::Float64)

Lineshape equation for Lorentzian (Cauchy) distribution.
"""
function lorentzian(energy::Float64, width::Float64)
    return width / (width^2 + energy^2)
end

"""
    lorentzian_integral(E_range::Tuple, N, W)

Analytic form for the indefinite integral for a lorentzian within the energy interval defined as (E_a, E_b).
"""
function lorentzian_integral(E_interval::Tuple, W)
    E_a, E_b = E_interval
    return 1/π * (atan(E_b/W) - atan(E_a/W))
end

"""
    lorentzianbath{T} <: BathFunction

Struct containing fields relating a quantum bath with density of states defined by a Lorentizian function.
"""
struct lorentzianbath{T} <: BathFunction
    bathfunction :: Vector{T}
    bathdegeneracy :: Vector{T}
    bathtype :: Symbol
    N :: Int64 # number of electronic states

end

"""
    lorentzianbath(discretisation::BathDiscretisationScheme, N::Int64; W=4.5)

Args:
- `discretisation`: chosed discretisation scheme for quantum bath
- `N`: number of electronic energy states described by the quantum bath
- `W`: width of lorenzian lineshape

Note:
Important distinction here between, `N`, the number of electronic energy states used here and `N` (or `M`) used elsewhere to describe the number of discrete states.
This is before the discretisation and concerns how the quantum bath is constructed. It is something we set manually here for a Lorentzian bath, but would be calculated from an actual DOS if using that approach.
When in the wide band limit and we assume no state degeneracy, number of electronic states is equivalent to the number of discrete states, hence these terms have been used interchangably elsewhere in the code.
"""
function lorentzianbath(discretisation::BathDiscretisationScheme, N::Int64; W=4.5)
    W_au = austrip(W*u"eV")
    bathfunction = lorentzian.(discretisation.bathstates, W_au)
    bathdegeneracy = N .* lorentzian_integral.(discretisation_energy_intervals(discretisation), W_au)
    bathtype = :lorentzianbath
    return lorentzianbath(bathfunction, bathdegeneracy, bathtype, N)
end