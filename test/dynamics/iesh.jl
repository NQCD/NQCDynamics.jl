using Test
using NonadiabaticMolecularDynamics
using LinearAlgebra
using Random
using Distributions

n_states = 20
temperature = 9.5e-4
model = MiaoSubotnik(M=n_states-1)
atoms = Atoms(2000)
r = zeros(1, 1)
v = randn(1,1) * sqrt(5*temperature / atoms.masses[1])
n_electrons = nstates(model) ÷ 2

sim = Simulation{IESH}(atoms, model; DoFs=1, n_electrons=n_electrons)
u = DynamicsVariables(sim, v, r)
sim.method.state .= u.state

@testset "compute_overlap!" begin
    S = zeros(n_electrons, n_electrons)
    state = collect(1:n_electrons)
    ψ = Dynamics.get_quantum_subsystem(u)
    Dynamics.compute_overlap!(S, ψ, state)
    @test S == I # Check overlap is identity

    # All electrons in state 1, so they all have unity overlap
    # with first electron wavefunction, which is in state 1.
    state = ones(Int, n_electrons)
    Dynamics.compute_overlap!(S, ψ, state)
    @test all(S[:,1] .== 1) # Check only first column ones
    @test all(S[:,2:end] .== 0)
end

@testset "calculate_Akj" begin
    S = zeros(n_electrons, n_electrons)
    state = collect(1:n_electrons)
    ψ = Dynamics.get_quantum_subsystem(u)
    Akj = Dynamics.calculate_Akj(S, ψ, 1.0, state)
    @test Akj ≈ 1.0
end

@testset "evaluate_hopping_probability!" begin
    Calculators.update_electronics!(sim.calculator, hcat(20.0))
    z = deepcopy(u)
    rand!(Dynamics.get_quantum_subsystem(z))
    Dynamics.evaluate_hopping_probability!(sim, z, 1.0)
end


BLAS.set_num_threads(1)

temp = 5*temperature
harm = Harmonic(m=model.m, ω=model.ω)
z = ComponentVector(v=zeros(1,1), r=hcat(10.0))
thermal = Simulation{Langevin}(atoms, harm; temperature=temp, DoFs=1, γ=0.1)
Δ = Dict(:X => 1.0)
sample = MetropolisHastings.run_monte_carlo_sampling(thermal, hcat(0.0), Δ, 1e6)

vel = Normal(0, sqrt(temp / atoms.masses[1]))

dist = InitialConditions.DynamicalDistribution(vel, sample.R, (1,1))

using Plots
plotly()

@time res = Ensembles.run_ensemble_standard_output(sim, (0.0, 5e4), dist;
    output=(:position, :population, :adiabatic_population), trajectories=2, abstol=1e-7, reltol=1e-4, saveat=0.0:50:5e4)

# This should plot the impurity population for each of the trajectories.
# It should look like the inset of fig 5.c in Miao Subotnik.
# Unfortunately it does not :(
impurity_population = zeros(length(0:50:5e4), length(res))
for (i, traj) in enumerate(res)
    # plot!(traj.t.*model.ω, 1 .- [p[1] for p in traj.population])
    impurity_population[:,i] .= 1 .- [p[1] for p in traj.population]
end

avg = mean(impurity_population, dims=2)
plot(res[1].t.*model.ω, avg)

    
