using NQCDynamics
using Test
using OrdinaryDiffEq: solve, MagnusGL4
using LinearAlgebra: norm

nstates = 100
tspan = (0.0, 100.0)
c0 = zeros(ComplexF64, nstates)
c0[1] = 1

@testset "Norm preservation" begin
    prob = DynamicsMethods.ElectronicODEProblem(c0, tspan, nstates)

    A = randn(nstates, nstates)
    A = (A * A')
    for i=1:nstates
        A[i,i] = 0
        for j=i+1:nstates
            A[j,i] *= -1
        end
    end

    DynamicsMethods.update_parameters!(prob.p, rand(nstates),
        [A for i=1:3, j=1:3],
        rand(3,3),
        tspan[2]
    )
    sol = solve(prob, MagnusGL4(krylov=true), dt=2.0)
    @test all(isapprox.(norm.(sol.u), 1))
end
