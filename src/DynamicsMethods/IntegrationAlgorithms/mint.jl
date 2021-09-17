
using UnPack: @unpack
using MuladdMacro: @muladd
using StaticArrays: SMatrix
using OrdinaryDiffEq: OrdinaryDiffEq
using LinearAlgebra: Hermitian, tr

using NonadiabaticMolecularDynamics.DynamicsMethods: MappingVariableMethods

"""
    MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm

Second order symplectic momentum integral algorithm applied to NRPMD.

# Reference

[J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
"""
struct MInt <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

OrdinaryDiffEq.isfsal(::MInt) = false

mutable struct MIntCache{uType,uEltypeNoUnits} <: OrdinaryDiffEq.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    cayley::Vector{Matrix{uEltypeNoUnits}}
end

function OrdinaryDiffEq.alg_cache(::MInt,u,rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
    tmp = zero(u)
    cayley = RingPolymers.cayley_propagator(p.beads, dt; half=true)
    MIntCache(u, uprev, tmp, cayley)
end

function OrdinaryDiffEq.initialize!(_, ::MIntCache) end

@muladd function OrdinaryDiffEq.perform_step!(integrator, cache::MIntCache, repeat_step=false)
    @unpack dt,uprev,u,p = integrator
    @unpack cayley = cache

    calc = p.calculator

    copyto!(u, uprev)
    
    RingPolymers.transform_to_normal_modes!(p.beads, u)
    step_C!(DynamicsUtils.get_velocities(u), DynamicsUtils.get_positions(u), cayley)
    RingPolymers.transform_from_normal_modes!(p.beads, u)

    Calculators.update_electronics!(calc, DynamicsUtils.get_positions(u))

    propagate_mapping_variables!(calc, u, dt)

    for i=1:length(p.beads)
        for j=1:length(p.atoms)
            for k=1:p.DoFs
                Γ = get_gamma(calc.adiabatic_derivative[k,j,i], calc.eigenvalues[i], dt)
                E = transform_matrix(calc.eigenvectors[i], Γ)

                Ξ = get_xi(calc.adiabatic_derivative[k,j,i], calc.eigenvalues[i], dt)
                F = transform_matrix(calc.eigenvectors[i], Ξ)

                qmap = MappingVariableMethods.get_mapping_positions(u, i)
                pmap = MappingVariableMethods.get_mapping_momenta(u, i)
                V = calc.derivative[k,j,i]
                force = get_mapping_nuclear_force(qmap, pmap, E, F, V, dt)
                DynamicsUtils.get_velocities(u)[k,j,i] -= force / p.atoms.masses[j]
            end
        end
    end


    RingPolymers.transform_to_normal_modes!(p.beads, u)
    step_C!(DynamicsUtils.get_velocities(u), DynamicsUtils.get_positions(u), cayley)
    RingPolymers.transform_from_normal_modes!(p.beads, u)

end

function step_C!(v::AbstractArray{T,3}, r::AbstractArray{T,3}, cayley::Vector{<:AbstractMatrix}) where {T}
    for I in CartesianIndices(v)
        rtmp = cayley[I[3]][1,1] * r[I] + cayley[I[3]][1,2] * v[I]
        vtmp = cayley[I[3]][2,1] * r[I] + cayley[I[3]][2,2] * v[I]
        r[I] = rtmp
        v[I] = vtmp
    end
end

function propagate_mapping_variables!(calc, u, dt)
    for i=1:length(calc.eigenvalues)
        C = get_C_propagator(calc.eigenvalues[i], calc.eigenvectors[i], dt)
        D = get_D_propagator(calc.eigenvalues[i], calc.eigenvectors[i], dt)
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
                           E::AbstractMatrix, F::AbstractMatrix,
                           diabatic_derivative::Hermitian, dt::Real)
    trVdt = tr(diabatic_derivative) * dt
    return 0.5 * (q'*E*q + p'*E*p - trVdt) - q'*F*p
end

DynamicsMethods.select_algorithm(::RingPolymerSimulation{<:MappingVariableMethods.NRPMD}) = MInt()
