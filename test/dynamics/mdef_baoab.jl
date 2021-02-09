using Test
using NonadiabaticMolecularDynamics
using StochasticDiffEq
using LinearAlgebra
using DiffEqNoiseProcess
using DiffEqDevTools
using Random
Random.seed!(100)

u0 = zeros(2,2)
v0 = ones(2,2)
γ = 1
f1_harmonic(v,u,p,t) = -u
f2_harmonic(v,u,p,t) = v
function g(u,p,t)
    Λ = diagm(ones(length(u)))
    σ = 0.2
    Λ, σ
end

ff_harmonic = DynamicalSDEFunction(f1_harmonic,f2_harmonic,g)
prob1 = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,0.5))

sol1 = solve(prob1,MDEF_BAOAB();dt=1/10,save_noise=true)

f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
function g_iip(du,u,p,t)
    Λ = diagm(ones(length(u)))
    σ = 0.2
    du.x[1] .= Λ
    du.x[2] .= σ
end

prob2 = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,0.5);
                           noise=NoiseWrapper(sol1.W))

sol2 = solve(prob2,MDEF_BAOAB();dt=1/10)
@test sol1[:] ≈ sol2[:]

dts = (1/2) .^ (8:-1:3)
res = analyticless_test_convergence(dts, prob1, MDEF_BAOAB(), (1/2)^12;
    trajectories=Int(1e2), use_noise_grid=false)
@test res.𝒪est[:weak_final] ≈ 1 atol=0.5
res = analyticless_test_convergence(dts, prob2, MDEF_BAOAB(), (1/2)^12;
    trajectories=Int(1e2), use_noise_grid=false)
@test res.𝒪est[:weak_final] ≈ 1 atol=0.5
