using NQCDynamics
using Test
using OrdinaryDiffEq: solve, MagnusGL4, Tsit5
using LinearAlgebra: norm, diagind

nstates = 10
tspan = (0.0, 100.0)
c0 = zeros(ComplexF64, nstates)
c0[1] = 1

A = randn(nstates, nstates)
A = (A * A')
for i=1:nstates
    A[i,i] = 0
    for j=i+1:nstates
        A[j,i] *= -1
    end
end

velocity = rand(3,3)
eigenvalues = rand(nstates)

@testset "ElectronicODEProblem, norm preservation" begin
    prob = DynamicsMethods.ElectronicODEProblem(c0, tspan, nstates)

    DynamicsMethods.update_parameters!(prob.p, eigenvalues,
        [A for i=1:3, j=1:3],
        velocity,
        tspan[2]
    )
    sol = solve(prob, MagnusGL4(krylov=true), dt=2.0)
    @test all(isapprox.(norm.(sol.u), 1))
end

@testset "DensityMatrixODEProblem, norm preservation" begin
    u0 = c0 * c0'
    prob = DynamicsMethods.DensityMatrixODEProblem(u0, tspan, nstates)

    DynamicsMethods.update_parameters!(prob.p, eigenvalues,
        [A for i=1:3, j=1:3],
        velocity,
        tspan[2]
    )
    sol = solve(prob, Tsit5())
    total_population = [real(sum(u[diagind(u)])) for u in sol.u]
    @test all(isapprox.(total_population, 1))
end
