using NonadiabaticMolecularDynamics
using Test
using ComponentArrays: ComponentVector

@testset "OutputFinal" begin
    output = Ensembles.OutputFinal()
    @test output(1:10, 1) == (10, false)
end

@testset "OutputDissociation" begin
    output = Ensembles.OutputDissociation(1.0, (1, 2))
    r = [1 0; 0 0; 0 0]
    v = zeros(3, 2)
    u = ComponentVector(v=v, r=r)
    sol = ComponentVector(u=[u])
    @test output(sol, 1) == (0, false)
    r = [2 0; 0 0; 0 0]
    u = ComponentVector(v=v, r=r)
    sol = ComponentVector(u=[u])
    @test output(sol, 1) == (1, false)
    @test Ensembles.output_template(output, u) == 0
end

@testset "OutputQuantisedDiatomic" begin
    sim = Simulation(Atoms(1), Harmonic())
    output = Ensembles.OutputQuantisedDiatomic(sim; height=1, normal_vector = [1, 1, 1])
    output = Ensembles.OutputQuantisedDiatomic(sim)
    u = DynamicsVariables(sim, hcat(1), hcat(1))
    @test Ensembles.output_template(output, u) == (0, 0)
end

@testset "PopulationCorrelationFunction" begin

    @testset "Ehrenfest" begin
        sim = Simulation{Ehrenfest}(Atoms(2000), TullyModelTwo())
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Adiabatic())
        u = DynamicsVariables(sim, hcat(20/2000), hcat(-5), SingleState(1, Adiabatic()))
        sol = run_trajectory(u, (0, 1000.0), sim)
        out, cont = output(sol, 1)
        @test out[1][1,1] ≈ 1
        @test out[1][2,2] ≈ 0
        @test [o[1,2] for o in out] ≈ [1 - o[1,1] for o in out]
    end

    @testset "FSSH" begin
        sim = Simulation{FSSH}(Atoms(2000), TullyModelTwo())
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Adiabatic())
        u = DynamicsVariables(sim, hcat(20/2000), hcat(-5), SingleState(1, Adiabatic()))
        sol = run_trajectory(u, (0, 1000.0), sim)
        out, cont = output(sol, 1)
        @test out[1][1,1] ≈ 1
        @test out[1][2,2] ≈ 0
        @test [o[1,2] for o in out] ≈ [1 - o[1,1] for o in out]
    end

    @testset "eCMM" begin
        sim = Simulation{eCMM}(Atoms(2000), TullyModelTwo())
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
        u = DynamicsVariables(sim, hcat(20/2000), hcat(-5), SingleState(1, Diabatic()))
        sol = run_trajectory(u, (0, 1000.0), sim)
        out, cont = output(sol, 1)
    end
    
end
