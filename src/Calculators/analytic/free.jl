export Free

struct Free <: AdiabaticCalculator

    n_states::UInt8
    potential::Function
    derivative::Function

    Free() = new(1, x->0, zero)
end
