using NonadiabaticMolecularDynamics
using Test
using FiniteDiff

function test_motion!(sim)
    f(x) = evaluate_hamiltonian(sim, x)

    v = randn(sim.DoFs, length(sim.atoms))
    r = randn(sim.DoFs, length(sim.atoms))
    u = ClassicalDynamicals(v, r)
    du = zero(u)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    Dynamics.motion!(du, u, sim, 0.0)

    # Rdot = dH/dP = dH 
    @test get_positions(du) ≈ grad.x.x[1]./sim.atoms.masses' rtol=1e-3
    # Vdot = Pdot / mass = -dH/dR / mass
    @test get_velocities(du) ≈ -grad.x.x[2]./sim.atoms.masses' rtol=1e-3
end

function test_velocity!(sim)
    u = randn(sim.DoFs, length(sim.atoms))
    v = randn(sim.DoFs, length(sim.atoms))
    du = zero(u)

    Dynamics.velocity!(du, v, u, sim, 0.0)

    @test du ≈ v
end

function test_acceleration!(sim)
    f(x) = evaluate_potential_energy(sim, x)

    r = randn(sim.DoFs, length(sim.atoms))
    v = randn(sim.DoFs, length(sim.atoms))
    dv = zero(v)

    grad = FiniteDiff.finite_difference_gradient(f, r)

    Dynamics.acceleration!(dv, v, r, sim, 0.0)

    @test dv ≈ -grad ./ sim.atoms.masses'
end

include("classical.jl")
# include("langevin.jl")
# include("mdef.jl")
include("fssh.jl")
include("fermionic.jl")
include("nrpmd.jl")

include("ensembles.jl")

include("saving_callbacks.jl")
include("cell_boundary_callback.jl")