
export GatesHollowayElbow

"""
    GatesHollowayElbow()

Simple two state elbow potential from Gates and Holloway:
[Journal of Electron Spectroscopy and Related Phenomena, 64/65 (1993) 633-639](https://doi.org/10.1016/0368-2048(93)80131-5)

Has two diabatic states each comprised of the sum of a Morse and a repulsive potential.
The coupling between them is an exponential function of `z` (distance from the surface).
"""
@with_kw struct GatesHollowayElbow <: DiabaticModel
    n_states::UInt8 = 2

    λ₁::Float64 = 3.5
    λ₂::Float64 = 3.5
    z₀::Float64 = 1.4
    x₀::Float64 = 0.6
    α::Float64 = 1.028
    d::Float64 = austrip(4.67u"eV")
    z12::Float64 = 0.5
    c::Float64 = austrip(0.5u"eV")
    γ::Float64 = 1.0
end

function potential!(model::GatesHollowayElbow, V, R)
    @unpack λ₁, λ₂, z₀, x₀, α, d, z12, c, γ = model

    repel(x, λ, d) = exp(-λ*(x+d))
    morse(x, d, α) = d*(1-exp(-α*x))^2

    x = R[1,1]
    z = R[1,2]

    V[1,1] = morse(x, d, α) + repel(z, λ₁, z₀)
    V[2,2] = morse(z, d, α) + repel(x, λ₂, x₀)
    V.data[1,2] = c * repel(z, γ,-z12)

    return V
end

function derivative!(model::GatesHollowayElbow, D, R)
    @unpack λ₁, λ₂, z₀, x₀, α, d, z12, c, γ = model

    drepel(x, λ, d) = -λ*exp(-λ*(x+d))
    dmorse(x, d, α) = 2α*d*(1-exp(-α*x))*exp(-α*x)

    x = R[1,1]
    z = R[1,2]

    # x derivative
    D[1,1][1,1] = dmorse(x, d, α)
    D[1,1][2,2] = drepel(x, λ₂, x₀)
    # D[1,1].data[1,2] = 0

    # z derivative
    D[1,2][1,1] = drepel(z, λ₁, z₀)
    D[1,2][2,2] = dmorse(z, d, α)
    D[1,2].data[1,2] = c * drepel(z, γ, -z12)

    return D
end
