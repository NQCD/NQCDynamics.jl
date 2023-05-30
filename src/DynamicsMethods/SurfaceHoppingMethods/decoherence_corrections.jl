
struct DecoherenceCorrectionNone end
apply_decoherence_correction!(args...) = nothing

"""
    DecoherenceCorrectionEDC{T}

Energy decoherence correction of Granucci and Persico in J. Chem. Phys. 126, 134114 (2007).
"""
struct DecoherenceCorrectionEDC{T}
    C::T
end

DecoherenceCorrectionEDC() = DecoherenceCorrectionEDC(0.1)

"""
Eq. 17 in J. Chem. Phys. 126, 134114 (2007)

Modify the wavefunction coefficients in `ψ` after a successful surface hop.
"""
function apply_decoherence_correction!(
    ψ::AbstractVector, correction::DecoherenceCorrectionEDC,
    occupied_state::Integer, dt::T, E::AbstractVector, Ekin::T
) where {T}

    C = correction.C
    unoccupied_state_norm = zero(T)
    for i in eachindex(ψ)
        if i != occupied_state
            τ = (1 + C/Ekin) / abs(E[i] - E[occupied_state])
            ψ[i] = ψ[i] * exp(-dt/τ)
            unoccupied_state_norm += abs2(ψ[i])
        end
    end

    Cₘ = ψ[occupied_state]
    ψ[occupied_state] = Cₘ * sqrt((1 - unoccupied_state_norm) / abs2(Cₘ))
    return nothing
end
