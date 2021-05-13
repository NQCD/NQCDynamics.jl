using FiniteDiff

get_blank(sim::Simulation) = randn(sim.DoFs, length(sim.atoms))
get_blank(sim::RingPolymerSimulation) = RingPolymerArray(randn(sim.DoFs, length(sim.atoms), length(sim.beads)))

function test_motion!(sim)
    f(x) = evaluate_hamiltonian(sim, x)

    v = get_blank(sim)
    r = get_blank(sim)
    u = ClassicalDynamicals(v, r)
    du = zero(u)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    Dynamics.motion!(du, u, sim, 0.0)

    # Rdot = dH/dP = dH/dV / mass
    @test get_positions(du) ≈ get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    # Vdot = Pdot / mass = -dH/dR / mass
    @test get_velocities(du) ≈ -get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
end

function test_velocity!(sim)
    u = get_blank(sim)
    v = get_blank(sim)
    du = zero(u)

    Dynamics.velocity!(du, v, u, sim, 0.0)

    @test du ≈ v
end

function test_acceleration!(sim)
    f(x) = evaluate_potential_energy(sim, x)

    r = get_blank(sim)
    v = get_blank(sim)
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    Dynamics.acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end

function test_acceleration!(sim::RingPolymerSimulation)
    f(x) = evaluate_potential_energy(sim, x)

    r = get_blank(sim)
    v = get_blank(sim)
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    Dynamics.ring_polymer_acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end
