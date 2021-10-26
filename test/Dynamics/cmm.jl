
using Test
using NonadiabaticMolecularDynamics
using FiniteDiff
using LinearAlgebra: norm

function test_motion!(sim::Simulation{<:eCMM}, u)
    f(x) = DynamicsUtils.classical_hamiltonian(sim, x)

    grad = FiniteDiff.finite_difference_gradient(f, u)

    du = zero(u)
    DynamicsMethods.motion!(du, u, sim, 0.0)

    @test DynamicsUtils.get_positions(du) ≈ DynamicsUtils.get_velocities(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsUtils.get_velocities(du) ≈ -DynamicsUtils.get_positions(grad) ./ sim.atoms.masses' rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_positions(du) ≈ DynamicsMethods.MappingVariableMethods.get_mapping_momenta(grad) rtol=1e-3
    @test DynamicsMethods.MappingVariableMethods.get_mapping_momenta(du) ≈ -DynamicsMethods.MappingVariableMethods.get_mapping_positions(grad) rtol=1e-3
end

sim = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=0.6)
sim1 = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=0.5)

v = randn(1,1)
r = randn(1,1)
u = DynamicsVariables(sim, v, r, SingleState(1))
u1 = DynamicsVariables(sim1, v, r, SingleState(1))

test_motion!(sim, u)
test_motion!(sim1, u1)

# sol = run_trajectory(u, (0, 100.0), sim; output=(:hamiltonian, :position, :population), reltol=1e-10, abstol=1e-10)
# @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-3
# @test all(isapprox.(sum.(sol.population), 1, rtol=1e-3))

# sol = run_trajectory(u1, (0, 100.0), sim1; output=(:hamiltonian, :position, :population), reltol=1e-10, abstol=1e-10)
# @test sol.hamiltonian[1] ≈ sol.hamiltonian[end] rtol=1e-3
# @test all(isapprox.(sum.(sol.population), 1, rtol=1e-3))

# @testset "Estimators.diabatic_population" begin
#     initial = Float64[]
#     final = Float64[]
#     n = 10
#     for i=1:n
#         u = DynamicsVariables(sim, 0, 0, 2)
#         push!(initial, DynamicsMethods.MappingVariableMethods.Kele(u.qmap, u.pmap, 1, sim.method.γ))
#         push!(final, DynamicsMethods.MappingVariableMethods.Keleinv(u.qmap, u.pmap, 1, sim.method.γ))
#     end
#     @show sum(initial) / n
#     @show sum(final) / n
#     total = initial .* final
#     @show sum(total) / n
#     @show Estimators.diabatic_population(sim, u)
# end

@testset "generate_random_points_on_nsphere" begin
    points = DynamicsMethods.MappingVariableMethods.generate_random_points_on_nsphere(10, 1)
    @test norm(points) ≈ 1
    points = DynamicsMethods.MappingVariableMethods.generate_random_points_on_nsphere(10, 10)
    @test norm(points) ≈ 10
end

function evaluate_correlation(γ, n)
    sim = Simulation{eCMM}(Atoms(1), DoubleWell(); γ=γ)
    initial = Float64[]
    final = Float64[]
    for i=1:n
        u = DynamicsVariables(sim, 0, 0, 2)
        push!(initial, DynamicsMethods.MappingVariableMethods.mapping_kernel(u.qmap, u.pmap, 1, sim.method.γ))
        push!(final, DynamicsMethods.MappingVariableMethods.inverse_mapping_kernel(u.qmap, u.pmap, 1, sim.method.γ))
    end
    total = initial .* final .* NonadiabaticModels.nstates(sim)
    return sum(total) / n
end

gams = -0.5:0.05:3.5
out = []
for γ in gams
    push!(out, evaluate_correlation(γ, 10000))
end
plot(gams, out)

# points = []
# for i=1:100
#     push!(points, DynamicsMethods.MappingVariableMethods.generate_random_points(2))
# end
# x = [p[1] for p in points] .* 5
# y = [p[2] for p in points] .* 5
# using Plots
# scatter(x, y)
