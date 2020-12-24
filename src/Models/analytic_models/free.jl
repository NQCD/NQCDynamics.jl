export Free

"""
A model that has zero external potential everywhere.
"""
struct Free <: AdiabaticModel

    n_states::UInt
    potential!::Function
    derivative!::Function
    
    clearit!(out, x) = fill!(out, 0.0)

    Free() = new(1, clearit!, clearit!)
end
