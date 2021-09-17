using UnPack: @unpack
using MuladdMacro: @muladd
using DiffEqBase: @..
using LinearAlgebra: Diagonal
using NonadiabaticMolecularDynamics: RingPolymers

export BCOCB

struct BCOCB <: StochasticDiffEq.StochasticDiffEqAlgorithm end

StochasticDiffEq.alg_compatible(prob::DiffEqBase.AbstractSDEProblem, alg::BCOCB) = true

struct BCOCBCache{U,A,K,uEltypeNoUnits,F} <: StochasticDiffEq.StochasticDiffEqMutableCache
    tmp::U
    rtmp::A
    vtmp::A
    k::K
    cayley::Vector{Matrix{uEltypeNoUnits}}
    half::uEltypeNoUnits
    friction::F
end

function StochasticDiffEq.alg_cache(::BCOCB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
    tmp = zero(u)
    rtmp = zero(u.x[1])
    vtmp = zero(rtmp)
    k = zeros(size(rtmp))

    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=true)
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

function FrictionCache(sim::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.MDEF}, dt) 
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

function FrictionCache(sim::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}, dt)
    γ0 = sim.method.γ
    γ = [γ0, 2sqrt.(2sim.beads.normal_mode_springs[2:end])...]
    c1 = exp.(-dt.*γ)
    c2 = sqrt.(1 .- c1.^2)
    sqrtmass = 1 ./ repeat(sqrt.(sim.atoms.masses'), sim.DoFs)
    σ = zero(sqrtmass)

    LangevinCache(c1, c2, sqrtmass, σ)
end

function StochasticDiffEq.initialize!(integrator, cache::BCOCBCache)
    @unpack t,uprev,p = integrator
    v = integrator.uprev.x[1]
    r = integrator.uprev.x[2]

    integrator.f.f1(cache.k,v,r,p,t)
end

@muladd function StochasticDiffEq.perform_step!(integrator,cache::BCOCBCache,f=integrator.f)
    @unpack t,dt,sqdt,uprev,u,p,W = integrator
    @unpack rtmp, vtmp, k, cayley, half, friction = cache

    v1 = uprev.x[1]
    r1 = uprev.x[2]

    step_B!(vtmp, v1, half*dt, k)
    @.. rtmp = r1

    RingPolymers.transform_to_normal_modes!(p.beads, rtmp)
    RingPolymers.transform_to_normal_modes!(p.beads, vtmp)

    step_C!(vtmp, rtmp, cayley)

    step_O!(friction, integrator, vtmp, rtmp, t+half*dt)

    step_C!(vtmp, rtmp, cayley)

    RingPolymers.transform_from_normal_modes!(p.beads, rtmp)
    RingPolymers.transform_from_normal_modes!(p.beads, vtmp)

    @.. u.x[2] = rtmp

    integrator.f.f1(k,vtmp,u.x[2],p,t+dt)
    step_B!(u.x[1], vtmp, half*dt, k)
end

function step_C!(v::R, r::R, cayley::Vector{<:Matrix}) where {R<:RingPolymers.RingPolymerArray}
    for i in axes(r, 3)
        for j in v.quantum_atoms
            for k in axes(r, 1)
                rtmp = cayley[i][1,1] * r[k,j,i] + cayley[i][1,2] * v[k,j,i]
                vtmp = cayley[i][2,1] * r[k,j,i] + cayley[i][2,2] * v[k,j,i]
                r[k,j,i] = rtmp
                v[k,j,i] = vtmp
            end
        end
    end

    for j in v.classical_atoms
        for k in axes(r, 1)
            rtmp = cayley[1][1,1] * r[k,j,1] + cayley[1][1,2] * v[k,j,1]
            vtmp = cayley[1][2,1] * r[k,j,1] + cayley[1][2,2] * v[k,j,1]
            r[k,j,1] = rtmp
            v[k,j,1] = vtmp
        end
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

function step_O!(friction::LangevinCache, integrator, v::R, r::R, t) where {R<:RingPolymers.RingPolymerArray}
    @unpack W, p, dt, sqdt = integrator
    @unpack c1, c2, sqrtmass, σ = friction

    @.. σ = sqrt(get_temperature(p, t)) * sqrtmass

    for i in axes(r, 3)
        for j in v.quantum_atoms
            @. v[:,j,i] = c1[i] * v[:,j,i] + c2[i] * σ[:,j] * W.dW[:,j,i] / sqdt
        end
    end

    for j in v.classical_atoms
        @. v[:,j,1] = c1[1] * v[:,j,1] + c2[1] * σ[:,j] * W.dW[:,j,1] / sqdt
    end
end

DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.AbstractMDEF}) = BCOCB()
DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:DynamicsMethods.ClassicalMethods.ThermalLangevin}) = BCOCB()
