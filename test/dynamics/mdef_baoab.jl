using Test
using StochasticDiffEq
using LinearAlgebra
using RecursiveArrayTools
using DiffEqNoiseProcess

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
prob = DynamicalSDEProblem(ff_harmonic,g,v0,u0,(0.0,5.0))

sol1 = solve(prob,MDEF_BAOAB();dt=1/10,save_noise=true)

f1_harmonic_iip(dv,v,u,p,t) = dv .= f1_harmonic(v,u,p,t)
f2_harmonic_iip(du,v,u,p,t) = du .= f2_harmonic(v,u,p,t)
function g_iip(du,u,p,t)
    Λ = diagm(ones(length(u)))
    σ = 0.2
    du.x[1] .= Λ
    du.x[2] .= σ
end

prob = DynamicalSDEProblem(f1_harmonic_iip,f2_harmonic_iip,g_iip,v0,u0,(0.0,5.0);
                           noise=NoiseWrapper(sol1.W))

sol2 = solve(prob,MDEF_BAOAB();dt=1/10)
@test sol1[:] ≈ sol2[:]