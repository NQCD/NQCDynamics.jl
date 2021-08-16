using Test
using NonadiabaticMolecularDynamics
using StatsBase: mean
using ComponentArrays
using Random
Random.seed!(1)

ats = [3, 1, 40]
Ds = [1, 3, 3]
Ts = [0.5, 1, 10]
beads = [1, 4, 2]

@testset "Classical" begin

    @testset "get_proposal" begin
        for (natoms, DoFs, T) in zip(ats, Ds, Ts)
            sim = Simulation(Atoms(rand(natoms)), Harmonic(); DoFs=DoFs, temperature=T)
            proposals = InitialConditions.get_proposal(sim, Dict(:X=>0.5), 1)
            @test length(proposals.proposal) == natoms * DoFs * 2
        end
    end

    @testset "Energy expectation" begin
        for (natoms, DoFs, T) in zip(ats, Ds, Ts)
            sim = Simulation(Atoms(rand(natoms)), Harmonic(); DoFs=DoFs, temperature=T)
            R0 = rand(DoFs, natoms)
            u0 = ComponentVector(v=zero(R0), r=R0)
            chain = InitialConditions.sample_configurations(sim, u0, 1e5, Dict(:X=>1); move_ratio=0.01)
            energy = evaluate_hamiltonian.(sim, chain)
            @test mean(energy) / (DoFs*natoms) ≈ T rtol=1e-1
        end
    end
end

@testset "Ring polymer" begin

    @testset "get_proposal" begin
        for (natoms, DoFs, T, nbeads) in zip(ats, Ds, Ts, beads)
            sim = RingPolymerSimulation(Atoms(vcat([:C], fill(:H, natoms-1))), Harmonic(), nbeads; DoFs=DoFs, temperature=T, quantum_nuclei=[:H])
            proposals = InitialConditions.get_proposal(sim, Dict(:H=>0.5, :C=>0.1), 1)
            @test length(proposals.proposal) == natoms * DoFs * nbeads * 2
        end
    end

    @testset "Energy expectation" begin
        for (natoms, DoFs, T, nbeads) in zip(ats, Ds, Ts, beads)
            sim = RingPolymerSimulation(Atoms(rand(natoms)), Harmonic(), nbeads; DoFs=DoFs, temperature=T)
            R0 = rand(DoFs, natoms, nbeads)
            u0 = ComponentVector(v=zero(R0), r=R0)
            chain = InitialConditions.sample_configurations(sim, u0, 1e5, Dict(:X=>10); move_ratio=0.01)
            potential = evaluate_hamiltonian.(sim, chain)
            @test mean(potential) / (DoFs*natoms*nbeads) ≈ nbeads*T rtol=1e-1
        end
    end
end
