
export DarlingHollowayElbow

"""
    DarlingHollowayElbow()

Adiabatic elbow potential from Darling and Holloway: [Faraday Discuss., 1993, 96, 43-54](https://doi.org/10.1039/FD9939600043)
"""
@with_kw struct DarlingHollowayElbow <: AdiabaticModel

    # Elbow parameters
    d::Float64 = austrip(4.76u"eV")
    α::Float64 = 1.028
    p::Float64 = 2.0
    zoff::Float64 = 0.3

    # Dissociation barrier parameters
    V₀::Float64 = austrip(0.7u"eV")
    βx::Float64 = 0.05
    xb::Float64 = 1.0
    βz::Float64 = 0.0
    zb::Float64 = 0.0

    # Physisorption parameters
    V₁::Float64 = austrip(0.85u"eV")
    βr::Float64 = 1.8
    zr::Float64 = 0.25
    Cvw::Float64 = austrip(5.306u"eV")
    kc::Float64 = 0.5
    zvw::Float64 = -1.9

    m₁::Float64 = 1.3
    m₂::Float64 = 0.0
    m₃::Float64 = -0.1

    w₁::Float64 = 1.5
    w₂::Float64 = 1.5
    w₃::Float64 = 2.0
end

function potential!(model::DarlingHollowayElbow, V, R)

    @unpack d, α, p, zoff, V₀, βx, xb, βz, zb, V₁, βr, zr, Cvw, zvw, kc, m₁, m₂, m₃, w₁, w₂, w₃ = model

    y(x,z) = (1 + 1/x^4)*(1 + 1/z^4)
    f(x,z) = (y(x,z) - 1)^(-1/4)
    s(a,b) = (1 + tanh(a*b)) / 2
    function xₑ(x,z)
        xp = x + p
        zp = z + p
        fp = f(xp,zp)
        return fp + f(xp-fp, zp-fp) - p
    end

    zf(x,z) = z - zoff * s(x - m₁, w₁)

    Velbow(x,z) = d*(1 - exp(-α*xₑ(x,z)))^2
    Vbarrier(x,z) = V₀*exp(-βx*(x-xb)^2 - βz*(zf(x,z)-zb)^2) * s(x-zf(x,z)-m₂,w₂)

    Vᵣ(z) = V₁*exp(-βr*(z-zr))
    fc(x) = 1 - (2x*(1+x)+1)*exp(-2x)
    Vvw(z) = - Cvw/(z - zvw)^3 * fc(kc*(z-zvw))
    Vphys(x,z) = (Vᵣ(zf(x,z)) + Vvw(zf(x,z))) * s(zf(x,z)-x-m₂, w₂) * s(zf(x,z)-m₃, w₃)

    total_V(x,z) = Velbow(x,z) + Vbarrier(x,z) + Vphys(x,z)

    x = R[1,1]
    Z = R[1,2]
    V[1] = total_V(x, Z)

    return V
end
