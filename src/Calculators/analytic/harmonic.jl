export Harmonic

struct Harmonic <: AdiabaticCalculator

    n_states::UInt8
    potential::Function
    derivative::Function
    
    potential(R) = sum(R.^2) / 2
    derivative(R) = R

    Harmonic() = new(1, potential, derivative)
end
