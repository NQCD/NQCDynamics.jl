
using Test
using NQCDynamics: FastDeterminant
using LinearAlgebra: det

n = 5

@testset "det! $T" for T in (Float32, Float64, ComplexF32, ComplexF64)
    A1 = rand(T, n, n)
    A2 = copy(A1)
    ipiv = zeros(Int, n)
    @test FastDeterminant.det!(A1, ipiv) ≈ det(A2)
end
