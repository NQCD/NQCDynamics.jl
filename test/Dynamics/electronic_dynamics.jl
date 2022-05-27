using NQCDynamics
using Test
using OrdinaryDiffEq: Tsit5, solve, LinearExponential, update_coefficients!
using LinearAlgebra: norm, diagind

nstates = 100
tspan = (0.0, 10.0)
c0 = zeros(ComplexF64, nstates)
c0[1] = 1

function prepare_parameters!(p::DynamicsMethods.SurfaceHoppingMethods.ElectronicParameters)
    DynamicsMethods.SurfaceHoppingMethods.set_eigenvalues!(p, rand(nstates))
    A = randn(nstates, nstates)
    A = (A * A')
    for i=1:nstates
        A[i,i] = 0
        for j=i+1:nstates
            A[j,i] *= -1
        end
    end
    DynamicsMethods.SurfaceHoppingMethods.set_dynamical_coupling!(p, A)
end

@testset "Norm preservation" begin
    prob = DynamicsMethods.SurfaceHoppingMethods.ElectronicODEProblem(c0, tspan, nstates)
    prepare_parameters!(prob.p)
    sol = solve(prob, LinearExponential(), dt=2.0)
    @test all(isapprox.(norm.(sol.u), 1))
end
