using Test
using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic
using RecursiveArrayTools
using LinearAlgebra: diag
using StatsBase
using StochasticDiffEq
using ComponentArrays
using NonadiabaticMolecularDynamics: DynamicsMethods, DynamicsUtils
using NonadiabaticMolecularDynamics.DynamicsMethods: ClassicalMethods, IntegrationAlgorithms

atoms = Atoms([:H, :C])
sim = RingPolymerSimulation{MDEF}(atoms, NonadiabaticModels.ConstantFriction(Free(), 1), 3; temperature=10u"K")

v = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), length(sim.beads)))
r = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), length(sim.beads)))

@testset "friction!" begin
    gtmp = zeros(length(r), length(r))
    gtmp = zeros(sim.DoFs*length(sim.atoms),sim.DoFs*length(sim.atoms),length(sim.beads))
    F = ClassicalMethods.friction!(gtmp, r, sim, 0.0)
    hmass = sim.calculator.model.γ/atoms.masses[1]
    cmass = sim.calculator.model.γ/atoms.masses[2]
    for i=1:length(sim.beads)
        @test diag(F[:,:,i]) ≈ [hmass, hmass, hmass, cmass, cmass, cmass]
    end
end

prob = DynamicsMethods.create_problem(ComponentVector(v=v,r=r), (0.0, 0.5), sim)

@testset "step_C!" begin
    dt = 0.5
    c = NonadiabaticMolecularDynamics.cayley_propagator(sim.beads, dt; half=true)
    ω_k = NonadiabaticMolecularDynamics.get_matsubara_frequencies(length(sim.beads), sim.beads.ω_n)
    vbefore = copy(v)
    rbefore = copy(r)
    IntegrationAlgorithms.step_C!(v,r,c)
    for I in CartesianIndices(r)
        A = [0 1; -ω_k[I[3]] 0]
        a = exp(A*dt/2) * [rbefore[I], vbefore[I]]
        @test a ≈ [r[I], v[I]] rtol=1e-4
    end
end

@testset "step_O!" begin
    cache = StochasticDiffEq.alg_cache(IntegrationAlgorithms.BCOCB(), prob, ArrayPartition(v,r),0,0,sim,0,0,0,Float64,0,0,0,0,0,0,Val{true})
    integrator = init(prob,IntegrationAlgorithms.BCOCB();dt=0.5)
    IntegrationAlgorithms.step_O!(cache.friction, integrator, v, r, 0.0)
end

@testset "ThermalLangevin" begin
    atoms = Atoms([:H, :C])
    sim = RingPolymerSimulation{ThermalLangevin}(atoms, NonadiabaticModels.Free(), 3; temperature=100u"K", γ=0.01, DoFs=1)

    v = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), length(sim.beads)))
    r = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), length(sim.beads)))

    sol = run_trajectory(ArrayPartition(v,r), (0.0, 1e4), sim; dt=1, output=(:kinetic))
    # @test mean(sol.kinetic) ≈ austrip(100u"K") * length(sim.beads) rtol=5e-1
end

@testset "MDEF" begin
    atoms = Atoms([:H, :C])
    sim = RingPolymerSimulation{MDEF}(atoms, RandomFriction(Harmonic()), 3; temperature=100u"K", DoFs=1)

    v = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), length(sim.beads)))
    r = RingPolymerArray(zeros(sim.DoFs, length(sim.atoms), length(sim.beads)))

    sol = run_trajectory(ArrayPartition(v,r), (0.0, 1e2), sim; dt=1, output=(:kinetic))
end
