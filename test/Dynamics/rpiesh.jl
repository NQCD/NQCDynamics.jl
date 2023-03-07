using Test
using NQCDynamics
using NQCDynamics: Calculators
using OrdinaryDiffEq: Tsit5, MagnusMidpoint
using Unitful, UnitfulAtomic

kT = 9.5e-4
M = 30 # number of bath states
Γ = 6.4e-3
W = 6Γ / 2 # bandwidth  parameter

basemodel = MiaoSubotnik(;Γ)
bath = TrapezoidalRule(M, -W, W)
model = AndersonHolstein(basemodel, bath)
atoms = Atoms(2000)
n_beads = 4
r = randn(1,1,n_beads)
v = randn(1,1,n_beads)
n_electrons = M ÷ 2

sim = RingPolymerSimulation{AdiabaticIESH}(atoms, model, n_beads)
u = DynamicsVariables(sim, v, r)
sim.method.state .= u.state
SurfaceHoppingMethods.set_unoccupied_states!(sim)

@testset "DynamicsVariables" begin
    model = AndersonHolstein(basemodel, bath; fermi_level=0.0u"eV")
    sim = RingPolymerSimulation{AdiabaticIESH}(atoms, model, n_beads)
    distribution = NQCDistributions.FermiDiracState(0.0u"eV", 0u"K")
    u = DynamicsVariables(sim, v, r, distribution)
    @test u.state == 1:n_electrons

    model = AndersonHolstein(basemodel, bath; fermi_level=0.1u"eV")
    sim = RingPolymerSimulation{AdiabaticIESH}(atoms, model, n_beads)
    distribution = NQCDistributions.FermiDiracState(0.0u"eV", 300u"K")
    @test_throws ErrorException DynamicsVariables(sim, v, r, distribution)

    model = AndersonHolstein(basemodel, bath; fermi_level=0.1u"eV")
    sim = RingPolymerSimulation{AdiabaticIESH}(atoms, model, n_beads)
    distribution = NQCDistributions.FermiDiracState(0.1u"eV", 300u"K")
    avg = zeros(nstates(model))
    samples = 5000
    for i=1:samples
        u = DynamicsVariables(sim, v, r, distribution)
        avg[u.state] .+= 1
    end
    avg ./= samples
    eigs = Calculators.get_centroid_eigen(sim.calculator, r)
    occupations = NQCDistributions.fermi.(eigs.values, distribution.fermi_level, distribution.β)       

    @test avg ≈ occupations atol=0.2
end

@testset "create_problem" begin
    sim = RingPolymerSimulation{AdiabaticIESH}(atoms, model, n_beads)
    DynamicsMethods.create_problem(u, (0.0, 1.0), sim)
    @test sim.method.state == 1:15
    @test sim.method.unoccupied == 16:31
end

@testset "algorithm comparison" begin
    sim = RingPolymerSimulation{AdiabaticIESH}(atoms, model, n_beads; temperature=300u"K")
    u = DynamicsVariables(sim, 10randn(1,1,n_beads) ./ atoms.masses[1], randn(1,1,n_beads))
    tspan = (0.0, 1000.0)
    dt = 1.0
    output = (OutputPosition, OutputVelocity, OutputQuantumSubsystem, OutputTotalEnergy)
    @time traj1 = run_dynamics(sim, tspan, u; dt, algorithm=DynamicsMethods.IntegrationAlgorithms.BCBWavefunction(), output)
    @time traj2 = run_dynamics(sim, tspan, u; algorithm=Tsit5(), saveat=tspan[1]:dt:tspan[2], output, reltol=1e-8, abstol=1e-8)

    display(plot(traj1, :OutputTotalEnergy))
    display(plot(traj2, :OutputTotalEnergy))
    # We cannot compare these when hopping is happening since the trajectories will be different. 
end
