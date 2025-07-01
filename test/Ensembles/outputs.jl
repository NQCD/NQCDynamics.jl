using NQCDynamics
using Test
using ComponentArrays: ComponentVector

@testset "OutputFinal" begin
    output = OutputFinal
    sol = ComponentVector(u=1:10)
    @test output(sol, 1) == 10
end

@testset "OutputDissociation" begin
    output = OutputDissociation(1.0, (1, 2))
    r = [1 0; 0 0; 0 0]
    v = zeros(3, 2)
    u = ComponentVector(v=v, r=r)
    sol = ComponentVector(u=[u])
    @test output(sol, 1) == 0
    r = [2 0; 0 0; 0 0]
    u = ComponentVector(v=v, r=r)
    sol = ComponentVector(u=[u])
    @test output(sol, 1) == 1
end

@testset "OutputQuantisedDiatomic" begin
    # Maybe tack this onto a dynamics unit test. 
    output = OutputQuantisedDiatomic(; height=1, normal_vector = [1, 1, 1])
    output = OutputQuantisedDiatomic()
end

@testset "PopulationCorrelationFunction" begin

    @testset "Ehrenfest" begin
        sim = Simulation{Ehrenfest}(Atoms(2000), TullyModelTwo())
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Adiabatic())
        u = DynamicsVariables(sim, hcat(20/2000), hcat(-5), PureState(1, Adiabatic()))
        sol = run_dynamics(sim, (0, 1000.0), u; output, dt=0.01)
        out = sol[:PopulationCorrelationFunction]
        @test out[1][1,1] ≈ 1
        @test out[1][2,2] ≈ 0
        @test [o[1,2] for o in out] ≈ [1 - o[1,1] for o in out]
    end

    @testset "FSSH" begin
        sim = Simulation{FSSH}(Atoms(2000), TullyModelTwo())
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Adiabatic())
        u = DynamicsVariables(sim, hcat(20/2000), hcat(-5), PureState(1, Adiabatic()))
        sol = run_dynamics(sim, (0, 1000.0), u; output)
        out = sol[:PopulationCorrelationFunction]
        @test out[1][1,1] ≈ 1
        @test out[1][2,2] ≈ 0
        @test [o[1,2] for o in out] ≈ [1 - o[1,1] for o in out]
    end

    @testset "eCMM" begin
        sim = Simulation{eCMM}(Atoms(2000), TullyModelTwo())
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
        u = DynamicsVariables(sim, hcat(20/2000), hcat(-5), PureState(1, Diabatic()))
        sol = run_dynamics(sim, (0, 1000.0), u; output)
        out = sol[:PopulationCorrelationFunction]
    end
    
end
