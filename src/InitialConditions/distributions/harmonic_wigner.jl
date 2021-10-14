
using Distributions: Normal

"""
    MomentumHarmonicWigner(ω, β)

Wigner distribution in a 1D harmonic potential for the momentum
"""
function MomentumHarmonicWigner(ω, β)
    σ = sqrt(Q(ω, β) / β)
    Normal(0, σ)
end

"""
    PositionHarmonicWigner(ω, β)

Wigner distribution in a 1D harmonic potential for the position
"""
function PositionHarmonicWigner(ω, β)
    σ = sqrt(Q(ω, β) / β) / ω  
    Normal(0, σ)
end

"""
    VelocityHarmonicWigner(ω, β, m)

Wigner distribution in a 1D harmonic potential for the velocity
"""
function VelocityHarmonicWigner(ω, β, m)
    σ = sqrt(Q(ω, β) / β) / m
    Normal(0, σ)
end

"Quantum corrector for the Wigner distribution"
function Q(ω, β)
    halfβħω = β*ω/2
    return halfβħω / tanh(halfβħω)
end
