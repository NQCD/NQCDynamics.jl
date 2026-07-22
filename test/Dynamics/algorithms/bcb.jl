
using Test
using NQCDynamics
using Unitful, UnitfulAtomic
using Distributions: Normal
using OrdinaryDiffEq: OrdinaryDiffEqCore, DynamicalODEProblem
using DiffEqDevTools: DiffEqDevTools
using RecursiveArrayTools: ArrayPartition

atoms = Atoms(:H)
model = Harmonic(m=atoms.masses[1])
n_beads = 10
sim = RingPolymerSimulation{Classical}(atoms, model, n_beads; temperature=300u"K")
v = VelocityBoltzmann(n_beads*300u"K", atoms.masses, size(sim)[1:2])
r = Normal(0.0, 0.01)
d = DynamicalDistribution(v, r, size(sim))
u0 = rand(d)

function analytic(u0,p,t)
    v0 = copy(DynamicsUtils.get_velocities(u0))
    r0 = copy(DynamicsUtils.get_positions(u0))

    RingPolymerArrays.transform_to_normal_modes!(v0, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(r0, p.beads.transformation)

    ω = sqrt.(RingPolymers.get_matsubara_frequencies(n_beads, p.beads.ω_n).^2 .+ model.ω^2)
    ω = reshape(ω, 1, 1, n_beads)
    r = @. r0 * cos(ω*t) + v0 / ω * sin(ω.*t)
    v = @. v0 * cos(ω*t) + -ω * r0 * sin(ω.*t)

    RingPolymerArrays.transform_from_normal_modes!(v, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(r, p.beads.transformation)
    return ArrayPartition(v, r)
end

tspan = (0.0, 10.0)
func = OrdinaryDiffEqCore.DynamicalODEFunction(DynamicsMethods.ClassicalMethods.acceleration!, DynamicsUtils.velocity!; analytic)
prob = DynamicalODEProblem(
    func, Array(DynamicsUtils.get_velocities(u0)), Array(DynamicsUtils.get_positions(u0)), tspan, sim,
)

dts = (1/2) .^ (4:-1:-1)
res = DiffEqDevTools.test_convergence(dts, prob, DynamicsMethods.IntegrationAlgorithms.BCB())
@test res.𝒪est[:l∞] ≈ 2 atol=0.4
