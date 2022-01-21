using NQCDynamics
using Test
using Unitful
using ComponentArrays

include("utils.jl")

atoms = Atoms([:H])
model = NQCModels.Harmonic()

@testset "Classical" begin
    sim = Simulation{Classical}(atoms, model)

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ComponentVector(v=v, r=r)

    sol = run_trajectory(u0, (0.0, 1000.0), sim; dt=0.1, output=(:hamiltonian))
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end

@testset "Ring polymer classical" begin
    sim = RingPolymerSimulation{Classical}(atoms, model, 10; temperature=10000u"K")

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ComponentVector(v=v, r=r)

    sol = run_trajectory(u0, (0.0, 1000.0), sim; dt=0.1, output=(:hamiltonian))
    @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-2
end
