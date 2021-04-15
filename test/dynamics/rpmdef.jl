using Test
using NonadiabaticMolecularDynamics
using Unitful
using UnitfulAtomic
using RecursiveArrayTools
using LinearAlgebra: diag
using StochasticDiffEq
using BenchmarkTools

atoms = Atoms([:H, :C])
sim = RingPolymerSimulation{MDEF}(atoms, Models.FreeConstantFriction(1), 3; temperature=10u"K")

v = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), length(sim.beads)))
r = RingPolymerArray(rand(sim.DoFs, length(sim.atoms), length(sim.beads)))

@testset "friction!" begin
    gtmp = zeros(length(r), length(r))
    gtmp = zeros(sim.DoFs*length(sim.atoms),sim.DoFs*length(sim.atoms),length(sim.beads))
    F = Dynamics.friction!(gtmp, r, sim, 0.0)
    hmass = sim.calculator.model.γ/atoms.masses[1]
    cmass = sim.calculator.model.γ/atoms.masses[2]
    for i=1:length(sim.beads)
        @test diag(F[:,:,i]) ≈ [hmass, hmass, hmass, cmass, cmass, cmass]
    end
end

prob = Dynamics.create_problem(ClassicalDynamicals(v,r), (0.0, 0.5), sim)

@testset "step_C!" begin
    dt = 0.5
    c = NonadiabaticMolecularDynamics.cayley_propagator(sim.beads, dt; half=true)
    ω_k = NonadiabaticMolecularDynamics.get_matsubara_frequencies(length(sim.beads), sim.beads.ω_n)
    vbefore = copy(v)
    rbefore = copy(r)
    Dynamics.step_C!(v,r,c)
    for I in CartesianIndices(r)
        A = [0 1; -ω_k[I[3]] 0]
        a = exp(A*dt/2) * [rbefore[I], vbefore[I]]
        @test a ≈ [r[I], v[I]] rtol=1e-4
    end
end

@testset "step_O!" begin
    cache = StochasticDiffEq.alg_cache(BCOCB(), prob, ArrayPartition(v,r),0,0,sim,0,0,0,Float64,0,0,0,0,0,0,Val{true})
    integrator = init(prob,BCOCB();dt=0.5)
    Dynamics.step_O!(cache.friction, integrator, v, r, 0.0)
end
