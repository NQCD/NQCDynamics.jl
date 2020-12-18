using Test
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.Dynamics.DynamicsMethods
using Unitful
using DifferentialEquations

@test MDEF{Float64}(10, 3) isa MDEF
atomic_parameters = Atoms.AtomicParameters(Atoms.InfiniteCell(), [:H, :C])
system = System{MDEF}(atomic_parameters, Models.Analytic.Free(), 10u"K")

R = zeros(n_DoF(system), n_atoms(system)) 
P = zeros(n_DoF(system), n_atoms(system)) 
u = Phasespace(R, P)
du = zero(u)

n = n_DoF(system)*n_atoms(system)*2
set_force!(du, u, system)

@testset "random_force!" begin
    # Test that only the bottom left of the matrix is filled
    blank = zeros(n, n)
    random_force!(blank, u, system, 1.0)
    @test all(blank[1:end, 1:n÷2] .== 0)
    @test all(blank[1:n÷2, 1:end] .== 0)
    @test all(blank[n÷2+1:end, n÷2+1:end] .!= 0)
end

prob = SDEProblem(differential!, random_force!, u, (0.0, 1.0), system; noise_rate_prototype=zeros(n, n))
sol = solve(prob)