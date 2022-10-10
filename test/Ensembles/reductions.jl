
using NQCDynamics
using Test
using OrdinaryDiffEq

@testset "$reduction" for reduction in (MeanReduction(), SumReduction())
    prob = ODEProblem((u,p,t) -> 1.01u, [0.5], (0.0,1.0))
    prob_func(prob, i, repeat) = remake(prob,u0=[rand()])
    output = Ensembles.EnsembleSaver((OutputFinal,))
    ensemble = EnsembleProblem(prob, prob_func=prob_func, output_func=output, reduction=reduction)
    sol = solve(ensemble, Tsit5(), trajectories=1e3, batch_size=20)
    if reduction isa Ensembles.SumReduction
        @test sol.u[:OutputFinal][1] / 1e3 ≈ 0.5 * exp(1.01) rtol=1e-1
    else
        @test sol.u[:OutputFinal][1] ≈ 0.5 * exp(1.01) rtol=1e-1
    end
end
