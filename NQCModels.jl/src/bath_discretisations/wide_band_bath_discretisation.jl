
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

function setcoupling!(out::AbstractVector, bathcoupling::AbstractVector, coupling::Real, couplings_rescale::Real=1.0)
    @inbounds for i in eachindex(out)
        out[i] = bathcoupling[i] * coupling * couplings_rescale
    end
end

function setcoupling!(out::AbstractVector, bathcoupling::Real, coupling::Real, couplings_rescale::Real=1.0)
    fill!(out, bathcoupling * coupling * couplings_rescale)
end
