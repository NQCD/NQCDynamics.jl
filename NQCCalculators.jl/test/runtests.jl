using ParallelTestRunner
using LinearAlgebra
using RingPolymerArrays
using NQCModels
using NQCCalculators

# Allocation-based tests aren't great - change these values if there's a good reason to.
# Eigenvalue allocations limit: Regardless of maths backend and Julia version, eigenvalues should take fewer allocations than this or something might be wrong
const MAX_EIGEN_ALLOCATIONS = 60000

runtests(NQCCalculators, ARGS)
