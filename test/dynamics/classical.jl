using NonadiabaticMolecularDynamics
using Test
using Unitful

include("utils.jl")

atoms = Atoms([:H])
model = Models.Harmonic()

@testset "Classical" begin
    sim = Simulation{Classical}(atoms, model; DoFs=1)

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ClassicalDynamicals(v, r)

    sol = Dynamics.run_trajectory(u0, (0.0, 1000.0), sim; dt=0.1, output=(:hamiltonian))
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end

@testset "Ring polymer classical" begin
    sim = RingPolymerSimulation{Classical}(atoms, model, 10; DoFs=1, temperature=10000u"K")

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ClassicalDynamicals(v, r)

    sol = Dynamics.run_trajectory(u0, (0.0, 1000.0), sim; dt=0.1, output=(:hamiltonian))
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end
