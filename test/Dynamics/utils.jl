using FiniteDiff
using ComponentArrays: ComponentVector
using NQCDynamics: DynamicsMethods, DynamicsUtils, nbeads, Estimators
using RingPolymerArrays: RingPolymerArray

get_blank(sim::Simulation) = randn(size(sim))
get_blank(sim::RingPolymerSimulation) = RingPolymerArray(randn(size(sim)))

function test_motion!(sim)
    f(x) = DynamicsUtils.classical_hamiltonian(sim, x)

    v = get_blank(sim)
    r = get_blank(sim)
    u = ComponentVector(v=v, r=r)
    du = zero(u)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    DynamicsMethods.motion!(du, u, sim, 0.0)

    # Rdot = dH/dP = dH/dV / mass
    @test DynamicsUtils.get_positions(du) ≈ DynamicsUtils.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    # Vdot = Pdot / mass = -dH/dR / mass
    @test DynamicsUtils.get_velocities(du) ≈ -DynamicsUtils.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
end

function test_velocity!(sim)
    r = get_blank(sim)
    v = get_blank(sim)
    du = zero(v)

    DynamicsUtils.velocity!(du, v, r, sim, 0.0)

    @test du ≈ v
end

function test_acceleration!(sim)
    f(x) = DynamicsUtils.classical_potential_energy(sim, x)

    r = get_blank(sim)
    v = get_blank(sim)
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    DynamicsMethods.ClassicalMethods.acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end

function test_acceleration!(sim::RingPolymerSimulation)
    f(x) = DynamicsUtils.classical_potential_energy(sim, x)

    r = get_blank(sim)
    v = get_blank(sim)
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    DynamicsMethods.ClassicalMethods.ring_polymer_acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end
