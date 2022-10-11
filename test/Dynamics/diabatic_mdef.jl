using Test
using NQCDynamics
using Unitful, UnitfulAtomic

atoms = Atoms(1u"u")
nstates = 100
x = austrip.((1.0:0.1:5)u"Å")
Γ = 0.25 # 0.02, 0.1, 0.25, 0.5, 1.0, 2.0
ΔE = Γ / 140
T = 5ΔE
thossmodel = ErpenbeckThoss(;Γ)
model = WideBandBath(thossmodel; step=ΔE, bandmin=-fld(nstates,2)*ΔE, bandmax=ΔE*fld(nstates,2))

@testset "Friction comparison" begin
    function fric(r, sim)
        r = [r;;]
        Λ = zeros(1, 1)
        t = 0.0
        ClassicalMethods.friction!(Λ, r, sim, t)
        return Λ[1]
    end

    σ = 20ΔE
    sim = Simulation{DiabaticMDEF}(atoms, model; temperature=T, friction_method=ClassicalMethods.DirectQuadrature(model.ρ, 1/T))
    friction_dq = fric.(x, sim)
    sim = Simulation{DiabaticMDEF}(atoms, model; temperature=T, friction_method=ClassicalMethods.WideBandExact(model.ρ, 1/T))
    friction_wb = fric.(x, sim)
    sim = Simulation{DiabaticMDEF}(atoms, model; temperature=T, friction_method=ClassicalMethods.GaussianBroadening(σ, 1/T))
    friction_gb = fric.(x, sim)
    sim = Simulation{DiabaticMDEF}(atoms, model; temperature=T, friction_method=ClassicalMethods.OffDiagonalGaussianBroadening(σ, 1/T))
    friction_ongb = fric.(x, sim)

    @test friction_dq ≈ friction_wb rtol=1e-1
    @test friction_gb ≈ friction_ongb rtol=1e-1
    @test friction_dq ≈ friction_ongb rtol=1e-1
end

@testset "Simulation{DiabaticMDEF}" begin
    sim = Simulation{DiabaticMDEF}(atoms, model; temperature=T, friction_method=ClassicalMethods.WideBandExact(model.ρ, 1/T))
    u = DynamicsVariables(sim, rand(1,1), rand(1,1))
    run_dynamics(sim, (0.0, 1.0), u; dt=0.1, output=OutputDynamicsVariables)
end

@testset "RingPolymerSimulation{DiabaticMDEF}" begin
    n_beads = 10
    sim = RingPolymerSimulation{DiabaticMDEF}(atoms, model, n_beads; temperature=T, friction_method=ClassicalMethods.WideBandExact(model.ρ, 1/T))
    u = DynamicsVariables(sim, rand(1,1,n_beads), rand(1,1,n_beads))
    run_dynamics(sim, (0.0, 1.0), u; dt=0.1, output=OutputDynamicsVariables)
end
