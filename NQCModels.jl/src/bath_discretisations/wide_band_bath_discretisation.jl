# ----------------------------------------- Will be made redundant ----------------------------------------- #
abstract type WideBandBathDiscretisation end
NQCModels.nstates(bath::WideBandBathDiscretisation) = length(bath.bathstates)

function fillbathstates!(out::Hermitian, bath::WideBandBathDiscretisation)
    diagonal = view(out, diagind(out)[2:end])
    copy!(diagonal, bath.bathstates)
end

function fillbathcoupling!(out::Hermitian, coupling::Real, bath::WideBandBathDiscretisation, couplings_rescale::Real=1.0)
    first_row = @view out.data[1, 2:end] 
    setcoupling!(first_row, bath.bathcoupling, coupling, couplings_rescale)

    return nothing
end

function fillbathstates!(out::Hermitian, bath::BathDiscretisationScheme)
    diagonal = view(out, diagind(out)[2:end])
    copy!(diagonal, bath.bathstates)
end

function fillbathcoupling!(out::Hermitian, coupling::Real, bath::BathDiscretisationScheme, couplings_rescale::Real=1.0)
    first_row = @view out.data[1, 2:end] 
    setcoupling!(first_row, bath.bathcoupling, coupling, couplings_rescale)

    return nothing
end

function setcoupling!(out::AbstractVector, bathcoupling::AbstractVector, coupling::Real, couplings_rescale::Real=1.0)
    @inbounds for i in eachindex(out)
        out[i] = bathcoupling[i] * coupling * couplings_rescale
    end
end

function setcoupling!(out::AbstractVector, bathcoupling::Real, coupling::Real, couplings_rescale::Real=1.0)
    fill!(out, bathcoupling * coupling * couplings_rescale)
end
# ---------------------------------------------------------------------------------------------------------- #

struct widebandbath{T} <: BathFunction
    bathfunction :: Vector{T}
    bathdegeneracy :: Vector{T}
    bathtype :: Symbol
    N :: Int64
end

"""
    widebandbath(discretisation::BathDiscretisationScheme)

Wide band bath function, bathdegeneracy set to unity - Density of States ignored.
Same as old baths.
"""
function widebandbath(discretisation::BathDiscretisationScheme)
    bathfunction = ones(NQCModels.nstates(discretisation))
    bathdegeneracy = bathfunction
    bathtype = :widebandbath
    return widebandbath(bathfunction, bathdegeneracy, bathtype, NQCModels.nstates(discretisation))
end

"""
    widebandbath(discretisation::BathDiscretisationScheme, N::Int64)

Wide band bath function, bath degeneracy is a constant scaled by `N/nstates(disretisation)`
"""
function widebandbath(discretisation::BathDiscretisationScheme, N::Int64)
    bathfunction = ones(NQCModels.nstates(discretisation))
    bathdegeneracy = bathfunction .* N/NQCModels.nstates(discretisation)
    bathtype = :widebandbath
    return widebandbath(bathfunction, bathdegeneracy, bathtype, N)
end