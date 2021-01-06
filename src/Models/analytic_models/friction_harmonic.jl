using Random

export FrictionHarmonic

"""
The 1-state harmonic model
"""
struct FrictionHarmonic <: FrictionModel

    n_states::UInt
    potential!::Function
    derivative!::Function
    friction!::Function

    function FrictionHarmonic(mass=1, omega=1, r₀=0)

        potential!(V, R) = V .= sum(0.5 * mass * omega ^ 2 .* (R .- r₀) .^2)
        derivative!(D, R) = D .= mass * omega ^ 2 .* (R .- r₀)
        friction!(F, R) = randn!(F)

        new(1, potential!, derivative!, friction!)
    end
end
