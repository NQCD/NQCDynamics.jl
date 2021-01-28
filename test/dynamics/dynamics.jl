using NonadiabaticMolecularDynamics
using Test
using FiniteDiff

function test_motion!(sim)
    function f(x)
        evaluate_hamiltonian(sim, x)
    end

    R = randn(sim.DoFs, length(sim.atoms))
    P = randn(sim.DoFs, length(sim.atoms))
    u = Phasespace(R, P)

    du = zero(u)
    Dynamics.motion!(du, u, sim, 0.0)
    grad = FiniteDiff.finite_difference_gradient(f, u)
    @test get_positions(du) ≈ get_momenta(grad) rtol=1e-3
    @test get_momenta(du) ≈ -get_positions(grad) rtol=1e-3
end

include("langevin.jl")
include("mdef.jl")
include("fssh.jl")
include("fermionic.jl")
include("ensembles.jl")
include("nrpmd.jl")

include("saving_callbacks.jl")