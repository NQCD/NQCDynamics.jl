using UnPack, MuladdMacro
using DiffEqBase: @..
using StochasticDiffEq
using StochasticDiffEq: StochasticDiffEqMutableCache
import StochasticDiffEq: alg_compatible, alg_cache, initialize!, perform_step!

export BCOCB

struct BCOCB <: StochasticDiffEqAlgorithm end

alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::BCOCB) = true

struct BCOCBCache{U,A,K,uEltypeNoUnits,F} <: StochasticDiffEqMutableCache
    tmp::U
    rtmp::A
    vtmp::A
    k::K
    cayley::Vector{Matrix{uEltypeNoUnits}}
    half::uEltypeNoUnits
    friction::F
end

function alg_cache(::BCOCB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
    tmp = zero(u)
    rtmp = zero(u.x[1])
    vtmp = zero(rtmp)
    k = zeros(size(rtmp))

    cayley = NonadiabaticMolecularDynamics.cayley_propagator(p.beads, dt; half=true)
    half = uEltypeNoUnits(1//2)

    friction = FrictionCache(p, dt)

    BCOCBCache(tmp, rtmp, vtmp, k, cayley, half, friction)
end

struct MDEFCache{flatV,G,N,C,M,S}
    flatvtmp::flatV
    tmp1::flatV
    tmp2::flatV
    gtmp::G
    noise::N
    c1::C
    c2::C
    sqrtmass::M
    σ::S
end

function FrictionCache(sim::RingPolymerSimulation{<:MDEF}, dt) 
    DoFs = sim.DoFs
    atoms = length(sim.atoms)
    beads = length(sim.beads)

    flatvtmp = zeros(DoFs*atoms)
    tmp1 = zero(flatvtmp)
    tmp2 = zero(flatvtmp)

    gtmp = zeros(DoFs*atoms, DoFs*atoms, beads)
    noise = zero(flatvtmp)

    c1 = Diagonal(zeros(DoFs*atoms, DoFs*atoms))
    c2 = zero(c1)
    
    sqrtmass = 1 ./ sqrt.(repeat(sim.atoms.masses; inner=DoFs))
    σ = zero(sqrtmass)

    MDEFCache(flatvtmp, tmp1, tmp2, gtmp, noise, c1, c2, sqrtmass, σ)
end

struct LangevinCache{C,M,S}
    c1::C
    c2::C
    sqrtmass::M
    σ::S
end

function FrictionCache(sim::RingPolymerSimulation{<:Langevin}, dt)
    γ0 = sim.method.γ
    γ = [γ0, 2get_matsubara_frequencies(length(sim.beads), sim.beads.ω_n)...]
    c1 = exp.(-dt.*γ)
    c2 = sqrt.(1 - c1.^2)
    sqrtmass = 1 ./ sqrt.(repeat(sim.atoms.masses; inner=DoFs))
    σ = zero(sqrtmass)

    LangevinCache(c1, c2, sqrtmass, σ)
end

function initialize!(integrator, cache::BCOCBCache)
    @unpack t,uprev,p = integrator
    v = integrator.uprev.x[1]
    r = integrator.uprev.x[2]

    integrator.f.f1(cache.k,v,r,p,t)
end

@muladd function perform_step!(integrator,cache::BCOCBCache,f=integrator.f)
    @unpack t,dt,sqdt,uprev,u,p,W = integrator
    @unpack rtmp, vtmp, k, cayley, half, friction = cache

    v1 = uprev.x[1]
    r1 = uprev.x[2]

    step_B!(vtmp, v1, half*dt, k)
    @.. rtmp = r1

    transform!(rtmp, p.beads.U)
    transform!(vtmp, p.beads.U)

    step_C!(vtmp, rtmp, cayley)

    step_O!(friction, integrator, vtmp, rtmp, t+half*dt)

    step_C!(vtmp, rtmp, cayley)

    transform!(rtmp, p.beads.U)
    transform!(vtmp, p.beads.U)

    @.. u.x[2] = rtmp

    integrator.f.f1(k,vtmp,u.x[2],p,t+dt)
    step_B!(u.x[1], vtmp, half*dt, k)
end

function step_C!(v::R, r::R, cayley::Vector{<:Matrix}) where {R<:RingPolymerArray}
    for I in CartesianIndices(v)
        rtmp = cayley[I[3]][1,1] * r[I] + cayley[I[3]][1,2] * v[I]
        vtmp = cayley[I[3]][2,1] * r[I] + cayley[I[3]][2,2] * v[I]
        r[I] = rtmp
        v[I] = vtmp
    end
end

function step_O!(friction::MDEFCache, integrator, v, r, t)
    @unpack W, p, dt, sqdt = integrator
    @unpack flatvtmp, tmp1, tmp2, gtmp, noise, c1, c2, sqrtmass, σ = friction

    integrator.g(gtmp,r,p,t)
    Λ = gtmp
    @.. σ = sqrt(get_temperature(p, t)) * sqrtmass

    @views for i in axes(r, 3)
        @.. noise = σ * W.dW[:,:,i][:] / sqdt

        γ, c = LAPACK.syev!('V', 'U', Λ[:,:,i]) # symmetric eigen
        clamp!(γ, 0, Inf)
        for (j,i) in enumerate(diagind(c1))
            c1[i] = exp(-γ[j]*dt)
            c2[i] = sqrt(1 - c1[i]^2)
        end

        copyto!(flatvtmp, v[:,:,i])
        mul!(tmp1, transpose(c), flatvtmp)
        mul!(tmp2, c1, tmp1)
        mul!(flatvtmp, c, tmp2)
        
        mul!(tmp1, transpose(c), noise)
        mul!(tmp2, c2, tmp1)
        mul!(tmp1, c, tmp2)

        @.. flatvtmp += tmp1
        copyto!(v[:,:,i], flatvtmp)
    end

end

function step_O!(friction::LangevinCache, integrator, v, r, t)
    @unpack W, p, dt, sqdt = integrator
    @unpack c1, c2, sqrtmass, σ = friction

    @.. σ = sqrt(get_temperature(p, t)) * sqrtmass

    @views for i in axes(r, 3)
        @.. v[:,:,i] = c1[i] * v[:,:,i] + c2[i] * σ * W.dW[:,:,i] / sqdt
    end

end
