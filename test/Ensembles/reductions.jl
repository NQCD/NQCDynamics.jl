
using NonadiabaticMolecularDynamics
using Test
using OrdinaryDiffEq

@testset "select_reduction" begin
    @test Ensembles.select_reduction(:mean) isa Ensembles.MeanReduction
    @test Ensembles.select_reduction(:sum) isa Ensembles.SumReduction
    @test Ensembles.select_reduction(:append) isa Function
end

@testset "get_u_init" begin

    @testset "AbstractOutput, reduction=$reduction" for reduction in (:mean, :sum)
        reduction = Ensembles.select_reduction(reduction)
        stripped_kwargs = Dict()
        tspan = (0.0, 1.0)
        u0 = 1.0
        output = Ensembles.OutputFinal()
        u_init = Ensembles.get_u_init(reduction, stripped_kwargs, tspan, u0, output)
        @test u_init ≈ zero(u0)
    end

    @testset "reduction=:append" begin
        reduction = Ensembles.select_reduction(:append)
        output = Ensembles.OutputFinal()
        u0 = 1
        u_init = Ensembles.get_u_init(reduction, Dict(), (0, 1), u0, output)
        @test u_init == []
    end

    @testset "TimeCorrelationFunction, reduction=$reduction" for reduction in (:mean, :sum)
        reduction = Ensembles.select_reduction(reduction)
        sim = Simulation{FSSH}(Atoms(1), DoubleWell())
        stripped_kwargs = Dict(:saveat=>0.1)
        tspan = (0.0, 1.0)
        u0 = DynamicsVariables(sim, 1, 0, SingleState(1, Adiabatic()))
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
        u_init = Ensembles.get_u_init(reduction, stripped_kwargs, tspan, u0, output)
        @test u_init ≈ [zeros(2, 2) for _ ∈ tspan[1]:0.1:tspan[2]]
    end

end

@testset "$reduction" for reduction in (Ensembles.MeanReduction(0), Ensembles.SumReduction())
    prob = ODEProblem((u,p,t) -> 1.01u, 0.5, (0.0,1.0))
    prob_func(prob, i, repeat) = remake(prob,u0=rand())
    output = Ensembles.OutputFinal()
    u_init = Ensembles.output_template(output, rand())
    ensemble = EnsembleProblem(prob, prob_func=prob_func, output_func=output, reduction=reduction, u_init=u_init)
    sol = solve(ensemble, Tsit5(), trajectories=1e3, batch_size=20)
    if reduction isa Ensembles.SumReduction
        @test sol.u / 1e3 ≈ 0.5 * exp(1.01) rtol=1e-1
    else
        @test sol.u ≈ 0.5 * exp(1.01) rtol=1e-1
    end
end
