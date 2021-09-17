using FiniteDiff
using ComponentArrays: ComponentVector

get_blank(sim::Simulation) = randn(sim.DoFs, length(sim.atoms))
get_blank(sim::RingPolymerSimulation) = RingPolymerArray(randn(sim.DoFs, length(sim.atoms), length(sim.beads)))

function test_motion!(sim)
    f(x) = NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim, x)

    v = get_blank(sim)
    r = get_blank(sim)
    u = ComponentVector(v=v, r=r)
    du = zero(u)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    Dynamics.motion!(du, u, sim, 0.0)

    # Rdot = dH/dP = dH/dV / mass
    @test Dynamics.get_positions(du) ≈ Dynamics.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    # Vdot = Pdot / mass = -dH/dR / mass
    @test Dynamics.get_velocities(du) ≈ -Dynamics.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
end

function test_velocity!(sim)
    r = get_blank(sim)
    v = get_blank(sim)
    du = zero(v)

    Dynamics.DynamicsUtils.velocity!(du, v, r, sim, 0.0)

    @test du ≈ v
end

function test_acceleration!(sim)
    f(x) = NonadiabaticMolecularDynamics.evaluate_potential_energy(sim, x)

    r = get_blank(sim)
    v = get_blank(sim)
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    Dynamics.ClassicalMethods.acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end

function test_acceleration!(sim::RingPolymerSimulation)
    f(x) = NonadiabaticMolecularDynamics.evaluate_potential_energy(sim, x)

    r = get_blank(sim)
    v = get_blank(sim)
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    Dynamics.ClassicalMethods.ring_polymer_acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end
