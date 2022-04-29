
using NQCDynamics
using Test
using OrdinaryDiffEq

@testset "Reduction" begin
    @test Ensembles.Reduction(:mean) isa Ensembles.MeanReduction
    @test Ensembles.Reduction(:sum) isa Ensembles.SumReduction
    @test Ensembles.Reduction(:append) isa Ensembles.AppendReduction
end

@testset "get_u_init" begin

    @testset "Output, reduction=$reduction" for reduction in (:mean, :sum)
        reduction = Ensembles.Reduction(reduction)
        tspan = (0.0, 1.0)
        u0 = 1.0
        saveat = 1
        output = Ensembles.OutputFinal()
        u_init = Ensembles.get_u_init(reduction, saveat, Dict(), tspan, u0, output)
        @test u_init ≈ zero(u0)
    end

    @testset "reduction=:append" begin
        reduction = Ensembles.Reduction(:append)
        output = Ensembles.OutputFinal()
        u0 = 1
        saveat = 1
        u_init = Ensembles.get_u_init(reduction, saveat, Dict(), (0, 1), u0, output)
        @test u_init === nothing
    end

    @testset "TimeCorrelationFunction, reduction=$reduction" for reduction in (:mean, :sum)
        reduction = Ensembles.Reduction(reduction)
        sim = Simulation{FSSH}(Atoms(1), DoubleWell())
        saveat = 0.1
        tspan = (0.0, 1.0)
        u0 = DynamicsVariables(sim, 1, 0, PureState(1, Adiabatic()))
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
        u_init = Ensembles.get_u_init(reduction, saveat, Dict(), tspan, u0, output)
        @test u_init ≈ [zeros(2, 2) for _ ∈ tspan[1]:0.1:tspan[2]]
    end

end

@testset "$reduction" for reduction in (Ensembles.MeanReduction(0), Ensembles.SumReduction())
    prob = ODEProblem((u,p,t) -> 1.01u, [0.5], (0.0,1.0))
    prob_func(prob, i, repeat) = remake(prob,u0=[rand()])
    output = Ensembles.OutputFinal()
    ensemble = EnsembleProblem(prob, prob_func=prob_func, output_func=output, reduction=reduction, u_init=[0.0])
    sol = solve(ensemble, Tsit5(), trajectories=1e3, batch_size=20)
    if reduction isa Ensembles.SumReduction
        @test sol.u[1] / 1e3 ≈ 0.5 * exp(1.01) rtol=1e-1
    else
        @test sol.u[1] ≈ 0.5 * exp(1.01) rtol=1e-1
    end
end
