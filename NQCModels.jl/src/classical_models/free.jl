
"""
    Free()

Zero external potential everywhere. Useful for modelling free particles.

```jldoctest model
julia> model, R = Free(3), rand(3, 10);

julia> potential(model, R)
0.0

julia> derivative(model, R)
3×10 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
struct Free <: ClassicalModel
    dofs::Int
end

Free() = Free(1)

NQCModels.ndofs(free::Free) = free.dofs

# GPU extensions: ext/NQCModelsCUDAExt.jl:13, ext/NQCModelsMetalExt.jl:13
NQCModels.potential(::Free, ::AbstractMatrix{T}) where {T} = zero(T)

function NQCModels.potential!(::Free, V::Matrix{<:Number}, ::AbstractMatrix{T}) where {T} 
    V .= zero(T)
end

# GPU extensions: ext/NQCModelsCUDAExt.jl:15, ext/NQCModelsMetalExt.jl:15
function NQCModels.derivative!(::Free, out::AbstractMatrix, ::AbstractMatrix{T}) where {T}
    out .= 0
    fill!(out, zero(T))
end