
using UnPack: @unpack
using MuladdMacro: @muladd
using StaticArrays: SMatrix
using LinearAlgebra: Hermitian, tr
using OrdinaryDiffEq.OrdinaryDiffEqCore: get_fsalfirstlast
using NQCDynamics.DynamicsMethods: MappingVariableMethods
using NQCModels: nstates

OrdinaryDiffEq.isfsal(::RingPolymerMInt) = false

mutable struct RingPolymerMIntCache{uType,uEltypeNoUnits} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    cayley::Vector{Matrix{uEltypeNoUnits}}
end

OrdinaryDiffEq.get_fsalfirstlast(cache::RingPolymerMIntCache, u::Any) = (nothing, nothing)

function OrdinaryDiffEq.alg_cache(::RingPolymerMInt,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=true)
    RingPolymerMIntCache(u, uprev, tmp, cayley)
end

function OrdinaryDiffEq.initialize!(_, ::RingPolymerMIntCache) end

@muladd function OrdinaryDiffEq.perform_step!(integrator, integrator_cache::RingPolymerMIntCache, repeat_step=false)
    @unpack dt,uprev,u,p = integrator
    @unpack cayley = integrator_cache

    cache = p.cache

    copyto!(u, uprev)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    
    RingPolymerArrays.transform_to_normal_modes!(r, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(v, p.beads.transformation)
    step_C!(v, r, cayley)
    RingPolymerArrays.transform_from_normal_modes!(r, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(v, p.beads.transformation)

    NQCCalculators.update_cache!(cache, r)
    NQCCalculators.update_V̄!(cache, r)
    NQCCalculators.update_D̄!(cache, r)
    NQCCalculators.update_traceless_adiabatic_derivative!(cache, r)

    propagate_mapping_variables!(cache, u, dt)

    for i=1:nbeads(p)
        traceless_eigs = cache.eigen[i].values .- cache.V̄[i]
        for j=1:natoms(p)
            for k=1:ndofs(p)

                Γ = get_gamma(cache.traceless_adiabatic_derivative[k,j,i], traceless_eigs, dt)
                E = transform_matrix(cache.eigen[i].vectors, Γ)

                Ξ = get_xi(cache.traceless_adiabatic_derivative[k,j,i], traceless_eigs, dt)
                F = transform_matrix(cache.eigen[i].vectors, Ξ)

                qmap = MappingVariableMethods.get_mapping_positions(u, i)
                pmap = MappingVariableMethods.get_mapping_momenta(u, i)
                force = get_mapping_nuclear_force(qmap, pmap, E, F)
                DynamicsUtils.get_velocities(u)[k,j,i] -= force / p.atoms.masses[j]
                DynamicsUtils.get_velocities(u)[k,j,i] -= cache.D̄[k,j,i] / p.atoms.masses[j] * dt
            end
        end
    end

    RingPolymerArrays.transform_to_normal_modes!(r, p.beads.transformation)
    RingPolymerArrays.transform_to_normal_modes!(v, p.beads.transformation)
    step_C!(v, r, cayley)
    RingPolymerArrays.transform_from_normal_modes!(r, p.beads.transformation)
    RingPolymerArrays.transform_from_normal_modes!(v, p.beads.transformation)

end

function propagate_mapping_variables!(cache, u, dt)
    for i=1:length(cache.eigen)
        V̄ = tr(cache.potential[i]) / nstates(cache.model)
        traceless_eigs = cache.eigen[i].values .- V̄
        C = get_C_propagator(traceless_eigs, cache.eigen[i].vectors, dt)
        D = get_D_propagator(traceless_eigs, cache.eigen[i].vectors, dt)
        qmap = MappingVariableMethods.get_mapping_positions(u, i)
        pmap = MappingVariableMethods.get_mapping_momenta(u, i)
        q = C * qmap - D * pmap
        p = C * pmap + D * qmap
        qmap .= q
        pmap .= p
    end
end

"Get the `C` propagator for the mapping variables."
function get_C_propagator(eigenvalues::AbstractVector, eigenvectors::AbstractMatrix, dt::Real)
    cosine = Diagonal(cos.(eigenvalues .* dt))
    return eigenvectors * cosine * eigenvectors'
end

"Get the `D` propagator for the mapping variables."
function get_D_propagator(eigenvalues::AbstractVector, eigenvectors::AbstractMatrix, dt::Real)
    sine = Diagonal(sin.(-eigenvalues .* dt))
    return eigenvectors * sine * eigenvectors'
end

transform_matrix(U, M) = U * M * U'

"Get the `Γ` variable used to calculate the nuclear propagators."
function get_gamma(W::AbstractMatrix, Λ::AbstractVector, dt::Real)
    n = length(Λ)
    Γ = SMatrix{n,n}(
    (i != j ? sin((Λ[i] - Λ[j])*dt) * W[i,j] / (Λ[i] - Λ[j]) : W[i,j]*dt for j=1:n, i=1:n)
    )
    return Γ
end

"Get the `Ξ` variable used to calculate the nuclear propagators."
function get_xi(W::AbstractMatrix, Λ::AbstractVector, dt::Real)
    n = length(Λ)
    Ξ = SMatrix{n,n}(
    (i != j ? (1 - cos((Λ[i] - Λ[j])*dt)) * W[i,j] / (Λ[i] - Λ[j]) : 0 for j=1:n, i=1:n)
    )
    return Ξ
end

"Get the force due to the mapping variables."
function get_mapping_nuclear_force(q::AbstractVector, p::AbstractVector,
                           E::AbstractMatrix, F::AbstractMatrix)
    return 0.5 * (q'*E*q + p'*E*p) - q'*F*p
end
