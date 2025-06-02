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

    sol = run_dynamics(sim, (0.0, 1000.0), u0; dt=0.1, output=(OutputTotalEnergy))
    #@test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], rtol=1e-2)
end

@testset "Ring polymer classical" begin
    sim = RingPolymerSimulation{Classical}(atoms, model, 10; temperature=10000u"K")

    test_velocity!(sim)
    test_acceleration!(sim)
    test_motion!(sim)

    v = get_blank(sim)
    r = get_blank(sim)
    u0 = ComponentVector(v=v, r=r)

    sol = run_dynamics(sim, (0.0, 1000.0), u0; dt=0.1, output=(OutputTotalEnergy))
    #@test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], rtol=1e-2)
end

#= @testset "Fermion model classical adiabatic dynamics" begin
    model = WideBandBath(ErpenbeckThoss(;Γ=0.1u"eV"), bandmin=-5u"eV", bandmax=5u"eV", step=0.2u"eV")
    sim = Simulation{Classical}(atoms, model; temperature=300u"K")
    v = rand(VelocityBoltzmann(300u"K", atoms.masses, size(sim)))
    r = [model.model.morse.x₀;;]
    u0 = DynamicsVariables(sim, v, r)
    sol = run_dynamics(sim, (0.0, 900.0u"fs"), u0; dt=1u"fs", output=(OutputTotalEnergy))
    @test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], atol=1e-2)
end
 =#
#= @testset "Fermion model ring polymer adiabatic dynamics" begin
    model = WideBandBath(ErpenbeckThoss(;Γ=0.1u"eV"), bandmin=-5u"eV", bandmax=5u"eV", step=0.2u"eV")
    n_beads = 10
    sim = RingPolymerSimulation{Classical}(atoms, model, n_beads; temperature=300u"K")
    v = VelocityBoltzmann(n_beads*300u"K", atoms.masses, size(sim)[1:2])
    r = model.model.morse.x₀
    d = DynamicalDistribution(v, r, size(sim))
    u0 = rand(d)
    sol = run_dynamics(sim, (0.0, 900.0u"fs"), u0; dt=1u"fs", output=(OutputTotalEnergy))
    @test isapprox(sol[:OutputTotalEnergy][1], sol[:OutputTotalEnergy][end], atol=1e-2)
end =#
