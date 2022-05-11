
"""
    FastDeterminant

When computing many determinants in a loop the allocation and GC of the temporary arrays
for the pivots and workspace can contribute a large portion of the total runtime.

The functions here are copied from `LinearAlgebra` but modified so the user can provide the temprorary
arrays directly.
"""
module FastDeterminant

using Base: require_one_based_indexing
using LinearAlgebra: chkstride1, LU, det
using LinearAlgebra.BLAS: @blasfunc, BlasInt, BlasFloat
using LinearAlgebra.LAPACK: LAPACK, chkargsok, liblapack, getrf!

"""
    det!(A::AbstractMatrix{T}, ipiv::AbstractVector{BlasInt}) where {T<:BlasFloat}

Same as `det` but the user must provide `ipiv` and the input matrix `A` is overwritten.
"""
function det!(A::AbstractMatrix{T}, ipiv::AbstractVector{BlasInt}) where {T<:BlasFloat}
    return det(my_lu!(A, ipiv))
end

function my_lu!(A::StridedMatrix{T}, ipiv::AbstractVector{BlasInt}) where {T<:BlasFloat}
    lpt = getrf!(A, ipiv)
    return LU{T,typeof(A)}(lpt[1], lpt[2], lpt[3])
end

for (getrf, elty) in ((:dgetrf_,:Float64), (:sgetrf_,:Float32), (:zgetrf_,:ComplexF64), (:cgetrf_,:ComplexF32))
    @eval function LAPACK.getrf!(A::AbstractMatrix{$elty}, ipiv::AbstractVector{BlasInt})
        require_one_based_indexing(A)
        chkstride1(A)
        chkstride1(ipiv)
        m, n = size(A)
        lda  = max(1,stride(A, 2))
        length(ipiv) ≥ min(m,n) || throw(DimensionMismatch())
        info = Ref{BlasInt}()
        ccall((@blasfunc($getrf), liblapack), Cvoid,
              (Ref{BlasInt}, Ref{BlasInt}, Ptr{$elty},
               Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}),
              m, n, A, lda, ipiv, info)
        chkargsok(info[])
        A, ipiv, info[] #Error code is stored in LU factorization type
    end
end

end # module
