
"""
    FastDeterminant

When computing many determinants in a loop the allocation and GC of the temporary arrays
for the pivots and workspace can contribute a large portion of the total runtime.

Using FastLapackInterface we can reduce the allocations and improve the runtime performance.
"""
module FastDeterminant

using FastLapackInterface: LUWs
using LinearAlgebra: LAPACK, det, LU

"""
    det!(A::AbstractMatrix, ws::LUWs)

Same as `det` but the user must provide the `LU` workspace from FastLapackInterface.
"""
function det!(A::AbstractMatrix, ws::LUWs)
    return det(LU(LAPACK.getrf!(ws, A)...))
end

end # module
