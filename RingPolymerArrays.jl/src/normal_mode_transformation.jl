
"""
    NormalModeTransformation{T}

Normal mode transformation for ring polymer systems.

Contains the transformation matrix along with a temporary vector to perform allocation-free
transformations.
"""
struct NormalModeTransformation{T}
    U::Matrix{T}
    tmp::Vector{T}
    function NormalModeTransformation{T}(n) where {T}
        # Real and imaginary parts of the discrete Fourier transform matrix. 
        U = [
            k > n ÷ 2 ?
            sin(2π * j * k / n) : cos(2π * j * k / n)
            for j in 0:n-1, k in 0:n-1
        ]

        for col in eachcol(U)
            normalize!(col) # Normalize each column
        end

        return new{T}(U, zeros(n))
    end
end

"""
    transform_to_normal_modes!(A::AbstractArray{T,3}, transform::NormalModeTransformation) where {T}

Transform all degrees of freedom into normal mode coordinates from the primitive coordinates.

Assumes the array is not already in normal mode coordinates.
"""
function transform_to_normal_modes!(A::AbstractArray{T,3}, transform::NormalModeTransformation) where {T}
    transform!(A, transpose(transform.U), transform.tmp)
end

"""
    transform_from_normal_modes!(A::AbstractArray{T,3}, transform::NormalModeTransformation) where {T}

Transform all degrees of freedom from the normal mode coordinates into the primitive coordinates. 

Assumes the array is already in normal mode coordinates.
"""
function transform_from_normal_modes!(A::AbstractArray{T,3}, transform::NormalModeTransformation) where {T}
    transform!(A, transform.U, transform.tmp)
end

function transform!(A::AbstractArray{T,3}, U::AbstractMatrix{T}, tmp::AbstractVector{T}) where {T}
    @views for i in axes(A,2)
        for j in axes(A, 1)
            mul!(tmp, U, A[j,i,:])
            copy!(A[j,i,:], tmp)
        end
    end
    return A
end
