export DebyeSpinBoson
export DebyeBosonBath

"""
    DebyeSpinBoson(Nᵇ; ϵ=0, Δ=1, η=0.09, ωᶜ=2.5)

Spin boson model with `Nᵇ` bosons with Debye spectral density:
`J(ω) = η * ω * ωᶜ / (ω^2 + ωᶜ^2)`

# References
[P. Huo and D. F. Coker, Mol. Phys. 110, 1035 (2012)](https://doi.org/10.1080/00268976.2012.684896)
[Rekik et al. J. Chem. Phys. 138, 144106 (2013)](http://aip.scitation.org/doi/10.1063/1.4799272)
"""
struct DebyeSpinBoson{E,D,N,W} <: DiabaticModel
    n_states::UInt8
    ϵ::E
    Δ::D
    η::N
    ωᶜ::W
    ωⱼ::Vector{W}
    cⱼ::Vector{W}
end

function DebyeSpinBoson(Nᵇ; ϵ=0, Δ=1, η=0.09, ωᶜ=2.5)

    ωᵐ = 10ωᶜ
    c(ω) = sqrt(2η * atan(ωᵐ / ωᶜ) / (π * Nᵇ)) * ω
    ω(j) = tan(j * atan(ωᵐ / ωᶜ) / Nᵇ) * ωᶜ

    ωⱼ = ω.(1:Nᵇ)
    cⱼ = c.(ωⱼ)

    DebyeSpinBoson(UInt8(2), ϵ, Δ, η, ωᶜ, ωⱼ, cⱼ)
end

function potential!(model::DebyeSpinBoson, V::Hermitian, R::AbstractMatrix)
    r = @view R[1,:]

    v0 = 0.0
    for i in eachindex(model.ωⱼ)
        v0 += model.ωⱼ[i]^2 * r[i]^2 / 2
    end
    V[1,1] = v0 + model.ϵ
    V[2,2] = v0 - model.ϵ

    for i in eachindex(model.cⱼ)
        V[1,1] += model.cⱼ[i] * r[i]
        V[2,2] -= model.cⱼ[i] * r[i]
    end

    V.data[1,2] = model.Δ
end

function derivative!(model::DebyeSpinBoson, D::AbstractMatrix{<:Hermitian}, R::AbstractMatrix)
    r = @view R[1,:]

    for i in eachindex(r)
        d0 = model.ωⱼ[i]^2 * r[i]
        D[1,i][1,1] = d0 + model.cⱼ[i]
        D[1,i][2,2] = d0 - model.cⱼ[i]
    end
end

"""
    DebyeBosonBath(Nᵇ; ωᶜ=2.5)

Bosonic bath with Debye spectral density.
Useful for sampling the bath uncoupled from the spin for spin-boson dynamics.
"""
struct DebyeBosonBath{W} <: AdiabaticModel
    n_states::UInt8
    ωᶜ::W
    ωⱼ::Vector{W}
end

function DebyeBosonBath(Nᵇ; ωᶜ=2.5)

    ωᵐ = 10ωᶜ
    ω(j) = tan(j * atan(ωᵐ / ωᶜ) / Nᵇ) * ωᶜ
    ωⱼ = ω.(1:Nᵇ)

    DebyeBosonBath(UInt8(2), ωᶜ, ωⱼ)
end

function potential!(model::DebyeBosonBath, V::Vector, R::AbstractMatrix)
    r = @view R[1,:]
    sum!(V, model.ωⱼ .^2 .* r .^2 ./ 2)
end

function derivative!(model::DebyeBosonBath, D::AbstractMatrix, R::AbstractMatrix)
    r = @view R[1,:]

    for i in eachindex(r)
        D[1,i] = model.ωⱼ[i]^2 * r[i]
    end
end
