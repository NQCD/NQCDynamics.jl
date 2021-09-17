using Test
using NonadiabaticMolecularDynamics
using StochasticDiffEq
using LinearAlgebra
using DiffEqNoiseProcess
using DiffEqDevTools
using Random
using Logging
Random.seed!(100)

u0 = zeros(2,2)
v0 = ones(2,2)
γ = 1
f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
g(u,p,t) = diagm(ones(length(u)))

p = Simulation(Atoms([:H,:H]), NonadiabaticModels.Free(); DoFs=2, temperature=100)

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
prob1 = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,0.5),p)

sol1 = solve(prob1,DynamicsMethods.IntegrationAlgorithms.MDEF_BAOAB();dt=1/10,save_noise=true)

f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
g_iip(du,u,p,t) = du .= g(u,p,t)

prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,0.5),p;
                           noise=NoiseWrapper(sol1.W))

sol2 = solve(prob2,DynamicsMethods.IntegrationAlgorithms.MDEF_BAOAB();dt=1/10)
@test sol1[:] ≈ sol2[:]

disable_logging(Logging.Info)

dts = (1/2) .^ (8:-1:3)
res = analyticless_test_convergence(dts, prob1, DynamicsMethods.IntegrationAlgorithms.MDEF_BAOAB(), (1/2)^12;
    trajectories=Int(1e2), use_noise_grid=false)
@test res.𝒪est[:weak_final] ≈ 1 atol=0.5
