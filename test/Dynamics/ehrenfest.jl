using Test
using NQCDynamics
using NQCDynamics: DynamicsMethods, Calculators
using NQCDynamics.DynamicsMethods.EhrenfestMethods
using OrdinaryDiffEq
using Statistics: var
using DataFrames, CSV
using Interpolations
import JSON

benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "Ehrenfest Tests")

@test Ehrenfest{Float64}(2) isa Ehrenfest
atoms = Atoms(:H)

@testset "Ehrenfest" begin
    sim = Simulation{Ehrenfest}(atoms, NQCModels.DoubleWell())

    r = zeros(size(sim))
    v = randn(size(sim))
    u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
    du = zero(u)

    @test DynamicsUtils.get_quantum_subsystem(u) ≈ Complex.([1 0; 0 0])

    DynamicsMethods.motion!(du, u, sim, 0.0)

    Calculators.evaluate_nonadiabatic_coupling!(sim.calculator, r)
    @test sim.calculator.nonadiabatic_coupling ≈ -sim.calculator.nonadiabatic_coupling'

    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5())

    σ = DynamicsUtils.get_quantum_subsystem(u)
    dσ = zero(σ)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)

    @testset "get_diabatic_population" begin
        population = Estimators.diabatic_population(sim, u)
        @test population ≈ [0.5, 0.5]
    end

    @testset "Algorithm comparison" begin
        atoms = Atoms(2000)
        sim = Simulation{Ehrenfest}(atoms, NQCModels.TullyModelTwo())
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        output = (OutputTotalEnergy, OutputKineticEnergy, OutputPotentialEnergy, OutputPosition, OutputVelocity, OutputQuantumSubsystem, OutputAdiabaticPopulation)

        dt = 2e-3
        u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
        dyn_test = @timed run_dynamics(sim, (0.0, 50.0), u; output, dt)
        traj1 = dyn_test.value
        @test isapprox(var(traj1[:OutputTotalEnergy]), 0; atol=1e-6)
        benchmark_results["TullyModelTwo"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

        u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
        dyn_test = @timed run_dynamics(sim, (0.0, 50.0), u; algorithm=Feagin12(), output, reltol=1e-15, abstol=1e-15, saveat=dt)
        traj2 = dyn_test.value
        @test isapprox(var(traj2[:OutputTotalEnergy]), 0; atol=1e-6)
        benchmark_results["TullyModelTwo Feagin12"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

        @test traj1[:OutputKineticEnergy] ≈ traj2[:OutputKineticEnergy] rtol = 1e-3
        @test traj1[:OutputPotentialEnergy] ≈ traj2[:OutputPotentialEnergy] rtol = 1e-1
        @test traj1[:OutputVelocity] ≈ traj2[:OutputVelocity] rtol = 1e-3
        @test traj1[:OutputPosition] ≈ traj2[:OutputPosition] rtol = 1e-3
        @test traj1[:OutputQuantumSubsystem][end] ≈ traj2[:OutputQuantumSubsystem][end] rtol = 1e-1
    end

    @testset "FermiDirac Ehrenfest algorithm comparison" begin
        kT = 9.5e-4
        M = 30 # number of bath states
        Γ = 6.4e-3
        W = 6Γ / 2 # bandwidth  parameter

        basemodel = MiaoSubotnik(; Γ)
        bath = TrapezoidalRule(M, -W, W)
        model = AndersonHolstein(basemodel, bath)
        atoms = Atoms(2000)

        sim = Simulation{Ehrenfest}(atoms, model)
        v = zeros(1, 1)
        r = hcat(21.0)
        u = DynamicsVariables(sim, v, r, FermiDiracState(0.0, 0.0))
        tspan = (0.0, 2000.0)
        dt = 10.0
        output = (OutputTotalEnergy, OutputKineticEnergy, OutputPotentialEnergy, OutputPosition, OutputVelocity, OutputQuantumSubsystem)

        dyn_test = @timed run_dynamics(sim, tspan, u; dt, output, algorithm=Vern9(), abstol=1e-15, reltol=1e-15, saveat=dt)
        traj1 = dyn_test.value
        @test isapprox(var(traj1[:OutputTotalEnergy]), 0; atol=1e-6)
        benchmark_results["FermiDirac Vern9"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

        u = DynamicsVariables(sim, v, r, FermiDiracState(0.0, 0.0))
        
        dyn_test = @timed run_dynamics(sim, tspan, u; dt, output) # default algorithm is fixed timestep
        traj2 = dyn_test.value
        @test isapprox(var(traj2[:OutputTotalEnergy]), 0; atol=1e-6)
        benchmark_results["FermiDirac fixed timestep"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)

        @test traj1[:OutputKineticEnergy] ≈ traj2[:OutputKineticEnergy] rtol = 1e-3
        @test traj1[:OutputPotentialEnergy] ≈ traj2[:OutputPotentialEnergy] rtol = 1e-3
        @test traj1[:OutputVelocity] ≈ traj2[:OutputVelocity] rtol = 1e-3
        @test traj1[:OutputPosition] ≈ traj2[:OutputPosition] rtol = 1e-3
        @test traj1[:OutputQuantumSubsystem] ≈ traj2[:OutputQuantumSubsystem] rtol = 1e-2
    end

    @testset "Spin boson population dynamics" begin

        N = 100
        atoms = Atoms(fill(1, N))
        β = 5
        T = 1 / β
        density = OhmicSpectralDensity(2.5, 0.09)
        model = SpinBoson(density, N, 0.0, 1.0)

        position = reshape([PositionHarmonicWigner(ω, β, 1) for ω in model.ωⱼ], 1, :)
        velocity = reshape([VelocityHarmonicWigner(ω, β, 1) for ω in model.ωⱼ], 1, :)
        distribution = DynamicalDistribution(velocity, position, (1, 100)) * PureState(1)

        sim = Simulation{Ehrenfest}(atoms, model)

        saveat = 0:0.1:20
        trajectories = 500
        output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
        
        dyn_test = @timed begin
            run_dynamics(sim, (0.0, 20.0), distribution;
            saveat, trajectories, output, reduction=MeanReduction(), dt=0.1)
        end
        ensemble = dyn_test.value
        result = [p[1, 1] - p[1, 2] for p in ensemble[:PopulationCorrelationFunction]]

        data = CSV.read(joinpath(@__DIR__, "reference_data", "gao_saller_jctc_2020_fig2b.csv"), DataFrame; header=false)
        itp = linear_interpolation(data[!, 1], data[!, 2]; extrapolation_bc=Line())
        for i in eachindex(result)
            t = ensemble[:Time][i]
            true_value = itp(t)
            @test isapprox(true_value, result[i]; atol=0.2)
        end
        benchmark_results["SpinBoson population"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
    end
end

@testset "Ehrenfest RPMD" begin
    atoms = Atoms(2000)
    sim = RingPolymerSimulation{Ehrenfest}(atoms, NQCModels.TullyModelTwo(), 10)
    v = fill(100 / 2000, 1, 1, 10)
    r = fill(-10.0, 1, 1, 10)
    u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
    dyn_test = @timed run_dynamics(sim, (0.0, 500.0), u, output=OutputTotalEnergy, dt=0.01)
    solution =  dyn_test.value
    @test isapprox(var(solution[:OutputTotalEnergy]), 0; atol=1e-6)
    benchmark_results["Ehrenfest RPMD"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @Info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @Info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/ehrenfest.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
