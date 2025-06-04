using Test
using NQCDynamics
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra: Eigen
using Random: seed!
using NQCDynamics: DynamicsMethods, DynamicsUtils
using NQCCalculators
using NQCDynamics.SurfaceHoppingMethods: SurfaceHoppingMethods
using NQCDynamics.DynamicsUtils: get_positions, get_velocities
using Unitful
import JSON
benchmark_dir = get(ENV, "BENCHMARK_OUTPUT_DIR", "tmp/nqcd_benchmark")
benchmark_results = Dict{String, Any}("title_for_plotting" => "FSSH Tests")

@test SurfaceHoppingMethods.FSSH{Float64}(2, :standard) isa SurfaceHoppingMethods.FSSH
atoms = Atoms(2)

@testset "FSSH" begin
    sim = Simulation{FSSH}(atoms, NQCModels.DoubleWell())
    sim.method.state = 1

    r = zeros(size(sim))
    v = randn(size(sim))
    u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
    du = zero(u)
    
    update_cache!(sim.cache, r) # Required for the cache to initialise all fields
    @test DynamicsUtils.get_quantum_subsystem(u) ≈ Complex.([1 0; 0 0])

    DynamicsMethods.motion!(du, u, sim, 0.0)

    coupling = NQCCalculators.get_nonadiabatic_coupling(sim.cache, r)
    @test coupling ≈ -coupling'

    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback = SurfaceHoppingMethods.HoppingCallback)
    DynamicsMethods.SurfaceHoppingMethods.evaluate_hopping_probability!(
        sim,
        u,
        get_proposed_dt(integrator),
    )

    σ = DynamicsUtils.get_quantum_subsystem(u)
    dσ = zero(σ)
    DynamicsUtils.set_quantum_derivative!(dσ, u, sim)

    @testset "get_diabatic_population" begin
        NQCDynamics.NQCCalculators.update_cache!(sim.cache, r)
        population = Estimators.diabatic_population(sim, u)
        @test all(isapprox.(population, 0.5))
    end

    @testset "select_new_state" begin
        sim.method.hopping_probability .= [0, 1.0]
        u.state = 1
        @test 2 == SurfaceHoppingMethods.select_new_state(sim, u)
        sim.method.hopping_probability .= [0, 1.0]
        u.state = 2
        @test 2 == SurfaceHoppingMethods.select_new_state(sim, u)
        sim.method.hopping_probability .= [1.0, 0.0]
        u.state = 2
        @test 1 == SurfaceHoppingMethods.select_new_state(sim, u)
    end

    @testset "rescale_velocity!" begin
        get_velocities(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
        integrator.p.method.new_state = 2
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == false # Check hop is rejected 

        get_velocities(integrator.u) .= 1e5 # Set momentum to big to force hop
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == true # Check hop is accepted
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.cache.eigen.values .= [0.9, -0.3]
        integrator.p.cache.eigen.vectors .= [1 0; 0.0 1]
        eigs = DynamicsUtils.get_hopping_eigenvalues(
            integrator.p,
            DynamicsUtils.get_positions(integrator.u),
        )
        @test 1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 1, 2)
        @test -1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)
    end

    @testset "execute_hop!" begin
        r = get_positions(integrator.u)
        r .= rand()
        sim = integrator.p
        calc = sim.cache
        NQCCalculators.update_cache!(calc, r)
        DynamicsUtils.get_velocities(integrator.u) .= 2 # Set high momentum to ensure successful hop
        integrator.u.state = 1
        integrator.p.method.new_state = 2
        eigs = DynamicsUtils.get_hopping_eigenvalues(
            integrator.p,
            DynamicsUtils.get_positions(integrator.u),
        )
        ΔE = SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)

        KE_initial = DynamicsUtils.classical_kinetic_energy(
            integrator.p,
            get_velocities(integrator.u),
        )
        H_initial = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
        @test first(integrator.u.state) == 1

        SurfaceHoppingMethods.execute_hop!(integrator)

        KE_final = DynamicsUtils.classical_kinetic_energy(
            integrator.p,
            get_velocities(integrator.u),
        )
        H_final = DynamicsUtils.classical_hamiltonian(integrator.p, integrator.u)
        ΔKE = KE_final - KE_initial

        @test first(integrator.u.state) == 2 # Check state has changed
        @test H_final ≈ H_initial
        @test ΔKE ≈ -ΔE rtol = 1e-3 # Test for energy conservation
    end

    @testset "run_dynamics" begin
        atoms = Atoms(2000)
        sim = Simulation{FSSH}(atoms, NQCModels.TullyModelTwo())
        v = hcat(100 / 2000)
        r = hcat(-10.0)
        u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
        dyn_test = @timed run_dynamics(sim, (0.0, 500.0), u, output=OutputTotalEnergy, reltol=1e-6)
        solution =  dyn_test.value
        @test solution[:OutputTotalEnergy][1] ≈ solution[:OutputTotalEnergy][end] rtol=1e-2
        benchmark_results["FSSH"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
    end
end

@testset "RPSH" begin
    sim = RingPolymerSimulation{FSSH}(
        atoms,
        NQCModels.DoubleWell(),
        5;
        temperature = 1000u"K",
    )

    r = RingPolymerArray(randn(size(sim)))
    v = RingPolymerArray(randn(size(sim)))
    u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
    sim.method.state = first(u.state)

    NQCDynamics.NQCCalculators.update_cache!(sim.cache, r)
    problem = ODEProblem(DynamicsMethods.motion!, u, (0.0, 1.0), sim)
    integrator = init(problem, Tsit5(), callback = SurfaceHoppingMethods.HoppingCallback)
    SurfaceHoppingMethods.evaluate_hopping_probability!(sim, u, get_proposed_dt(integrator))

    @testset "rescale_velocity!" begin
        get_velocities(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
        integrator.p.method.new_state = 2
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == false # Check hop is rejected 

        get_velocities(integrator.u) .= 1e5 # Set momentum to big to force hop
        cont = SurfaceHoppingMethods.rescale_velocity!(integrator.p, integrator.u)
        @test cont == true # Check hop is accepted
    end

    @testset "calculate_potential_energy_change" begin
        integrator.p.cache.centroid_eigen.values .= [0.9, -0.3]
        integrator.p.cache.centroid_eigen.vectors .= [1 0; 0.0 1]
        eigs = DynamicsUtils.get_hopping_eigenvalues(
            integrator.p,
            DynamicsUtils.get_positions(integrator.u),
        )
        @test 1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 1, 2)
        @test -1.2 ≈ SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)
    end

    @testset "execute_hop!" begin
        NQCDynamics.NQCCalculators.update_cache!(integrator.p.cache, get_positions(integrator.u))
        get_velocities(integrator.u) .= 5 # Set high momentum to ensure successful hop
        integrator.u.state = [1] # Modifying SurfaceHoppingVariables, so supply a vector
        integrator.p.method.new_state = 2 # Modifying RPSH, so supply an Int. 
        eigs = DynamicsUtils.get_hopping_eigenvalues(
            integrator.p,
            DynamicsUtils.get_positions(integrator.u),
        )
        ΔE = SurfaceHoppingMethods.calculate_potential_energy_change(eigs, 2, 1)

        KE_initial =
            DynamicsUtils.centroid_classical_kinetic_energy(integrator.p, integrator.u)
        V_initial =
            DynamicsUtils.centroid_classical_potential_energy(integrator.p, integrator.u)
        H_initial = DynamicsUtils.centroid_classical_hamiltonian(integrator.p, integrator.u)

        @test first(integrator.u.state) == 1
        SurfaceHoppingMethods.execute_hop!(integrator)
        @test first(integrator.u.state) == 2 # Check state has changed

        KE_final =
            DynamicsUtils.centroid_classical_kinetic_energy(integrator.p, integrator.u)
        V_final =
            DynamicsUtils.centroid_classical_potential_energy(integrator.p, integrator.u)
        H_final = DynamicsUtils.centroid_classical_hamiltonian(integrator.p, integrator.u)

        @test H_initial ≈ H_final # Check ring polymer hamiltonian conserved
        ΔKE = (KE_final - KE_initial)
        ΔV = (V_final - V_initial)
        @test ΔKE ≈ -ΔE # Test for energy conservation
        @test ΔV ≈ -ΔKE
    end

    @testset "run_dynamics" begin
        atoms = Atoms(2000)
        sim = RingPolymerSimulation{FSSH}(atoms, TullyModelTwo(), 5; temperature = 0.01)
        v = fill(100 / 2000, size(sim))
        r = fill(-10.0, size(sim)) .+ randn(1, 1, 5)
        u = DynamicsVariables(sim, v, r, PureState(1, Adiabatic()))
        seed!(1)
        dyn_test = @timed run_dynamics(sim, (0.0, 1000.0), u, output=OutputTotalEnergy, dt=0.1)
        solution = dyn_test.value
        seed!(1)
        solution2 = run_dynamics(
            sim,
            (0.0, 1000.0),
            u,
            output = OutputTotalEnergy,
            algorithm = Tsit5(),
            abstol = 1e-10,
            reltol = 1e-10,
            saveat = 0:0.1:1000.0,
        )
        # Ring polymer Hamiltonian is not strictly conserved during hoppping
        @test solution[:OutputTotalEnergy][1] ≈ solution[:OutputTotalEnergy][end] rtol=1e-2
        benchmark_results["RPSH"] = Dict("Time" => dyn_test.time, "Allocs" => dyn_test.bytes)
    end

end

# Make benchmark directory if it doesn't already exist.
if !isdir(benchmark_dir)
    mkpath(benchmark_dir)
    @info "Benchmark data ouput directory created at $(benchmark_dir)."
else
    @info "Benchmark data ouput directory exists at $(benchmark_dir)."
end

# Output benchmarking dict
output_file = open("$(benchmark_dir)/FSSH.json", "w")
JSON.print(output_file, benchmark_results)
close(output_file)
