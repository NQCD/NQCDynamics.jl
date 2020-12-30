export Harmonic

"""
The 1-state harmonic model
"""
struct Harmonic <: AdiabaticModel

    n_states::UInt
    potential!::Function
    derivative!::Function

    function Harmonic(mass=1, omega=1, r₀=0)

        potential!(V, R) = V .= sum(0.5 * mass * omega ^ 2 .* (R .- r₀) .^2)
        derivative!(D, R) = D .= mass * omega ^ 2 .* (R .- r₀)

        new(1, potential!, derivative!)
    end
end
