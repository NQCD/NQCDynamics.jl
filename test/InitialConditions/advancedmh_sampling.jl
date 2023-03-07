using Test
using NQCDynamics
using NQCDynamics.InitialConditions
using StatsBase: mean
using ComponentArrays
using Random
using FiniteDiff
Random.seed!(1)

ats = [3, 1, 40]
Ds = [1, 3, 3]
Ts = [0.5, 1, 10]
beads = [1, 16, 2]

@testset "Classical" begin

    @testset "get_proposal" begin
        for (natoms, DoFs, T) in zip(ats, Ds, Ts)
            sim = Simulation(Atoms(rand(natoms)), Harmonic(dofs=DoFs); temperature=T)
            proposals = ThermalMonteCarlo.get_proposal(sim, Dict(:X=>0.5), 1, 0)
            @test length(proposals.proposal) == natoms * DoFs
        end
    end

    @testset "Energy expectation" begin
        for (natoms, DoFs, T) in zip(ats, Ds, Ts)
            sim = Simulation(Atoms(rand(natoms)), Harmonic(dofs=DoFs); temperature=T)
            R0 = randn(size(sim))
            chain = ThermalMonteCarlo.run_advancedmh_sampling(sim, R0, 1e4, Dict(:X=>1); move_ratio=0.5)
            energy = Estimators.@estimate potential_energy(sim, chain)
            @test energy / (DoFs*natoms) ≈ T/2 rtol=1e-1
        end
    end

    @testset "HMC gradient" begin
        for (natoms, DoFs, T) in zip(ats, Ds, Ts)
            sim = Simulation(Atoms(rand(natoms)), Harmonic(dofs=DoFs); temperature=T)
            density = ThermalMonteCarlo.get_density_function(sim)
            r = randn(prod(size(sim)))
            grad = FiniteDiff.finite_difference_gradient(density, r)
            grad_func = ThermalMonteCarlo.get_deriv_function(sim)
            @test grad ≈ grad_func(r)[2]
        end
    end

    @testset "Energy expectation HMC" begin
        for (natoms, DoFs, T) in zip(ats, Ds, Ts)
            sim = Simulation(Atoms(rand(natoms)), Harmonic(dofs=DoFs); temperature=T)
            R0 = randn(size(sim))
            chain, stats = ThermalMonteCarlo.run_advancedhmc_sampling(sim, R0, 1e4; verbose=false)
            energy = Estimators.@estimate potential_energy(sim, chain)
            @test energy / (DoFs*natoms) ≈ T/2 rtol=1e-1
        end
    end

end

@testset "Ring polymer" begin

    @testset "get_proposal" begin
        for (natoms, DoFs, T, nbeads) in zip(ats, Ds, Ts, beads)
            sim = RingPolymerSimulation(Atoms(vcat([:C], fill(:H, natoms-1))), Harmonic(dofs=DoFs), nbeads; temperature=T, quantum_nuclei=[:H])
            proposals = ThermalMonteCarlo.get_proposal(sim, Dict(:H=>0.5, :C=>0.1), 1, 0)
            @test length(proposals.proposal) == natoms * DoFs * nbeads
        end
    end

    @testset "Energy expectation" begin
        for (natoms, DoFs, T, nbeads) in zip(ats, Ds, Ts, beads)
            sim = RingPolymerSimulation(Atoms(rand(natoms)), Harmonic(dofs=DoFs), nbeads; temperature=T)
            R0 = zeros(DoFs, natoms, nbeads)
            chain = ThermalMonteCarlo.run_advancedmh_sampling(sim, R0, 1e4, Dict(:X=>1); move_ratio=0.3, internal_ratio=0.95)
            energy = DynamicsUtils.classical_potential_energy.(sim, chain) .+ DynamicsUtils.classical_spring_energy.(sim, chain)
            @test mean(energy) / (DoFs*natoms*nbeads) ≈ nbeads*T/2 rtol=1e-1
        end
    end

    @testset "Quantum harmonic oscillator" begin
        ω = 100
        T = 1
        sim = RingPolymerSimulation(Atoms(1), Harmonic(ω=ω), 200; temperature=T)
        r = zeros(size(sim))
        chain = ThermalMonteCarlo.run_advancedmh_sampling(sim, r, 1e5, Dict(:X=>1/100); move_ratio=0.0, internal_ratio=0.95)

        exact = ω/2 * coth(ω / 2T)
        estimate = Estimators.@estimate total_energy(sim, chain)
        @test exact ≈ estimate rtol=1e-1
    end

    @testset "Free ring polymer radius_of_gyration" begin
        T = 3.0
        analytic_Rg = sqrt(1/12T)
        sim = RingPolymerSimulation(Atoms(1), Free(), 10; temperature=T)
        r = zeros(size(sim))
        chain = ThermalMonteCarlo.run_advancedmh_sampling(sim, r, 4e4, Dict(:X=>1); move_ratio=0.0, internal_ratio=0.5)
        estimate = Estimators.@estimate radius_of_gyration(sim, chain)
        @test analytic_Rg ≈ estimate[1] atol=1e-1
    end

end
