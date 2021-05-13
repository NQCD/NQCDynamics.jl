export Free

"""
A model that has zero external potential everywhere.
"""
struct Free <: AdiabaticModel end

potential!(::Free, out::AbstractVector, ::AbstractMatrix) = out .= 0
derivative!(::Free, out::AbstractMatrix, ::AbstractMatrix) = out .= 0