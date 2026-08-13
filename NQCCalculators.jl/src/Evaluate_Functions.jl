#Definitions of evaluate functions

function evaluate_potential(cache::Abstract_ClassicalModel_Cache, r::AbstractMatrix)
    return NQCModels.potential(cache.model, r)
end

function evaluate_potential(cache::Abstract_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    potential = zero(cache.potential)
    @views @inbounds for i in beads(cache)
        potential[:,:,i] .= NQCModels.potential(cache.model, r[:,:,i])
    end
    return potential
end

function evaluate_potential(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    potential = zero(cache.potential)
    NQCModels.potential!(cache.model, potential, r)
    return potential
end

function evaluate_potential(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    potential = zero(cache.potential)
    @views @inbounds for i in beads(cache)
        NQCModels.potential!(cache.model, potential[i], r[:,:,i])
    end
    return potential
end

function evaluate_derivative(cache::Abstract_ClassicalModel_Cache, r::AbstractMatrix)
    derivative = zero(cache.derivative)
    NQCModels.derivative!(cache.model, derivative, r::AbstractMatrix)
    return derivative
end

function evaluate_derivative(cache::Abstract_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    derivative = zero(cache.derivative)
    @views @inbounds for i in beads(cache)
        NQCModels.derivative!(cache.model, derivative[:,:,i], r[:,:,i])
    end
    return derivative
end

function evaluate_derivative(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    derivative = zero(cache.derivative)
    NQCModels.derivative!(cache.model, derivative, r::AbstractMatrix)
    return derivative
end

function evaluate_derivative(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    derivative = zero(cache.derivative)
    @views @inbounds for i in beads(cache)
        NQCModels.derivative!(cache.model, derivative[:,:,i], r[:,:,i])
    end
    return derivative
end

function evaluate_friction(cache::Abstract_ClassicalModel_Cache, R::AbstractMatrix)
    friction = zero(cache.friction)
    NQCModels.friction!(cache.model, friction, R::AbstractMatrix)
    return friction
end

function evaluate_friction(cache::Abstract_ClassicalModel_Cache, R::AbstractArray{T,3}) where {T}
    friction = zero(cache.friction)
    @views for i in axes(R, 3)
        NQCModels.friction!(cache.model, friction[:,:,i], R[:,:,i])
    end
    return friction
end

function evaluate_friction(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    μ = NQCModels.fermilevel(cache.model)
    friction = zero(cache.friction)
    if cache.friction_method isa WideBandExact
        potential = evaluate_potential(cache, r)
        derivative = evaluate_derivative(cache, r)

        fill_friction_tensor!(friction, cache.friction_method, potential, derivative, r, μ)

        return friction
    elseif cache.friction_method isa Nothing
        return friction
    else
        ∂H = evaluate_adiabatic_derivative(cache, r)
        eigen = evaluate_eigen(cache, r)

        fill_friction_tensor!(friction, cache.friction_method, eigen, ∂H, r, μ)

        return friction
    end
end

function evaluate_friction(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    μ = NQCModels.fermilevel(cache.model)
    friction = zero(cache.friction)
    if cache.friction_method isa WideBandExact
        potential = evaluate_potential(cache, r)
        derivative = evaluate_derivative(cache, r)

        @views for i in axes(r, 3)
            fill_friction_tensor!(friction[:,:,i], cache.friction_method, potential[i], derivative[:,:,i], r[:,:,i], μ)
        end
        return friction
    elseif cache.friction_method isa Nothing
        return friction
    else
        ∂H = evaluate_adiabatic_derivative(cache, r)
        eigen = evaluate_eigen(cache, r)

        @views for i in axes(r, 3)
            fill_friction_tensor!(friction[:,:,i], cache.friction_method, eigen[i], ∂H[:,:,i], r[:,:,i], μ)
        end
        return friction
    end
end

function evaluate_centroid_friction(cache::Abstract_ClassicalModel_Cache, R::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(R)
    centroid_friction = zero(cache.centroid_friction)
    NQCModels.friction!(cache.model, centroid_friction, centroid)
    return centroid_friction
end

function evaluate_centroid_friction!(cache::Abstract_QuantumModel_Cache, R::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(R)
    μ = NQCModels.fermilevel(cache.model)
    centroid_friction = zero(cache.centroid_friction)

    if cache.friction_method isa WideBandExact
        potential = evaluate_centroid_potential(cache, R)
        derivative = evaluate_centroid_derivative(cache, R)

        fill_friction_tensor!(centroid_friction, cache.friction_method, potential, derivative, centroid, μ)

        return centroid_friction
    elseif cache.friction_method isa Nothing
        return centroid_friction
    else
        ∂H = evaluate_centroid_adiabatic_derivative(cache, R)
        eigen = evaluate_centroid_eigen(cache, R)

        fill_friction_tensor!(centroid_friction, cache.friction_method, ∂H, eigen, centroid, μ)

        return centroid_friction
    end
end

function correct_phase!(cache, eig::FastLapackInterface.HermitianEigenWs)
    @views for i in 1:length(eig.w)
        eig.Z[:,i] .*= 2*(0.5 - signbit(cache.phase_ref' * eig.Z[:,i]))
    end
    return nothing
end

function correct_phase!(cache, eig::AbstractMatrix)
    @views for i in 1:size(eig, 2)
        eig[:,i] .*= 2*(0.5 - signbit(cache.phase_ref' * eig[:,i]))
    end
    return nothing
end

function correct_phase!(phase_ref::AbstractVector, eig::AbstractMatrix)
    @views for i in 1:size(eig, 2)
        eig[:,i] .*= 2*(0.5 - signbit(phase_ref' * eig[:,i]))
    end
    return nothing
end

function correct_phase!(phase_ref::AbstractVector, eig::FastLapackInterface.HermitianEigenWs)
    @views for i in 1:length(eig.w)
        eig.Z[:,i] .*= 2*(0.5 - signbit(phase_ref' * eig.Z[:,i]))
    end
    return nothing
end

function evaluate_eigen(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    potential = evaluate_potential(cache, r)
    eig = HermitianEigenWs(zero(cache.eigen.Z) + I)
    FastLapackInterface.syevr!(eig, 'V', 'A', 'U', Matrix(potential) , 0.0, 0.0, 0, 0, 1e-12)
    correct_phase!(cache, eig)
    return eig
end

function evaluate_eigen(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    RP_eigen = [HermitianEigenWs(zero(cache.eigen[1].Z) + I, vecs=true) for _=1:length(beads(cache))]
    potential = evaluate_potential(cache, r)
    tmp = zeros(size(potential[1]))
    @inbounds for i in beads(cache)
        tmp .= potential[i]
        FastLapackInterface.syevr!(RP_eigen[i], 'V', 'A', 'U', tmp , 0.0, 0.0, 0, 0, 1e-12)
        correct_phase!(cache.phase_ref[i], RP_eigen[i])
    end
    return RP_eigen
end

function evaluate_adiabatic_derivative(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    diabatic_derivative = evaluate_derivative(cache, r)
    U = evaluate_eigen(cache, r)
    adiabatic_derivative = zero(cache.adiabatic_derivative)

    for I in eachindex(diabatic_derivative)
        adiabatic_derivative[I] .= U.Z' * diabatic_derivative[I] * U.Z
    end
    return adiabatic_derivative
end

function evaluate_adiabatic_derivative(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    derivative = evaluate_derivative(cache, r)
    eigen = evaluate_eigen(cache, r)
    adiabatic_derivative = zero(cache.adiabatic_derivative)

    for i in axes(derivative, 3) # Beads
        for j in axes(derivative, 2) # Atoms
            for k in axes(derivative, 1) # DoFs
                adiabatic_derivative[k,j,i] .= eigen[i].Z' * derivative[k,j,i] * eigen[i].Z
            end
        end
    end
    return adiabatic_derivative
end

function evaluate_centroid_adiabatic_derivative(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T, 3}) where {T}
    centroid_derivative = evaluate_centroid_derivative(cache, r)
    centroid_eigen = evaluate_centroid_eigen(cache, r)
    centroid_adiabatic_derivative = zero(cache.centroid_adiabatic_derivative)

    for I in eachindex(centroid_derivative)
        centroid_adiabatic_derivative[I] .= centroid_eigen.Z' * centroid_derivative[I] * centroid_eigen.Z
    end
    return centroid_adiabatic_derivative
end

function evaluate_inverse_difference_matrix!(out, eigenvalues)
    @inbounds for i in eachindex(eigenvalues)
        for j in eachindex(eigenvalues)
            out[j,i] = 1 / (eigenvalues[i] - eigenvalues[j])
        end
        out[i,i] = zero(eltype(out))
    end
    return nothing
end

"""
# References

- HammesSchifferTully_JChemPhys_101_4657_1994 Eq. (32)
- SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, section 2.3
"""
function evaluate_nonadiabatic_coupling(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    eigen = evaluate_eigen(cache, r)
    adiabatic_derivative = evaluate_adiabatic_derivative(cache, r)
    nonadiabatic_coupling = zero(cache.nonadiabatic_coupling)

    evaluate_inverse_difference_matrix!(cache.tmp_mat, eigen.w)

    @inbounds for I in NQCModels.dofs(cache)
        @. cache.nonadiabatic_coupling[I] = adiabatic_derivative[I] * cache.tmp_mat
    end
    return nonadiabatic_coupling
end

function evaluate_nonadiabatic_coupling(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    eigen = evaluate_eigen(cache, r)
    adiabatic_derivative = evaluate_adiabatic_derivative(cache, r)
    nonadiabatic_coupling = zero(cache.nonadiabatic_coupling)

    @inbounds for i in beads(cache)
        evaluate_inverse_difference_matrix!(cache.tmp_mat, eigen[i].w)
        for j in mobileatoms(cache)
            for k in dofs(cache)
                @. nonadiabatic_coupling[k,j,i] = adiabatic_derivative[k,j,i] * cache.tmp_mat
            end
        end
    end
    return nonadiabatic_coupling
end

function evaluate_centroid_nonadiabatic_coupling(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    eigen = evaluate_centroid_eigen(cache, r)
    centroid_adiabatic_derivative = evaluate_centroid_adiabatic_derivative(cache, r)
    centroid_nonadiabatic_coupling = zero(cache.centroid_nonadiabatic_coupling)

    evaluate_inverse_difference_matrix!(cache.tmp_mat, eigen.w)

    @inbounds for j in mobileatoms(cache)
        for k in dofs(cache)
            @. centroid_nonadiabatic_coupling[j,k] = centroid_adiabatic_derivative[j,k] * cache.tmp_mat
        end
    end
    return centroid_adiabatic_derivative
end





#RingPolymer specific functions


function evaluate_centroid(cache::Abstract_Cache, r::AbstractArray{T,3}) where {T}
    return RingPolymerArrays.get_centroid(r)
end

function evaluate_centroid_potential(cache::Abstract_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(r)
    return NQCModels.potential(cache.model, centroid)
end

function evaluate_centroid_potential(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(r)
    centroid_potential = zero(cache.centroid_potential)
    NQCModels.potential!(cache.model, centroid_potential, centroid)
    return centroid_potential
end

function evaluate_centroid_derivative(cache::Abstract_Cache, r::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(r)
    centroid_derivative = zero(cache.centroid_derivative)
    NQCModels.derivative!(cache.model, centroid_derivative, centroid)
    return centroid_derivative
end

function evaluate_centroid_eigen(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    potential = evaluate_centroid_potential(cache, r)
    tmp = zeros(size(potential))
    tmp .= potential
    eig = HermitianEigenWs(zero(cache.centroid_eigen.Z) + I)
    FastLapackInterface.syevr!(eig, 'V', 'A', 'U', tmp, 0.0, 0.0, 0, 0, 1e-12)
    correct_phase!(cache.centroid_phase_ref, eig)
    return eig
end

#evaluate tilde quantities

function evaluate_V̄(cache::Union{RingPolymer_QuantumModel_Cache, RingPolymer_QuantumFrictionModel_Cache}, r)
    potential = evaluate_potential(cache, r)
    V̄ = zero(cache.V̄)

    for i in 1:length(V̄)
        V̄[i] = tr(potential[i]) / nstates(cache.model)
    end
    return V̄
end

function evaluate_D̄(cache::Union{RingPolymer_QuantumModel_Cache,RingPolymer_QuantumFrictionModel_Cache}, r)
    derivative = evaluate_derivative(cache, r)
    D̄ = zero(cache.D̄)

    for I in eachindex(derivative)
        D̄[I] = tr(derivative[I]) / nstates(cache.model)
    end
    return D̄
end

#evaluate traceless quantities

function evaluate_traceless_potential(cache::Union{RingPolymer_QuantumModel_Cache, RingPolymer_QuantumFrictionModel_Cache}, r::AbstractArray{T,3}) where {T}
    potential = evaluate_potential(cache, r)
    V̄ = evaluate_V̄(cache, r)
    traceless_potential = zero(cache.traceless_potential)

    n = nstates(cache)
    for i in eachindex(potential)
        traceless_potential[i] = potential[i] - V̄[i].*I(n)
    end
    return traceless_potential
end

function evaluate_traceless_derivative(cache::Union{RingPolymer_QuantumModel_Cache, RingPolymer_QuantumFrictionModel_Cache}, r::AbstractArray{T,3}) where {T}
    derivative = evaluate_derivative(cache, r)
    D̄ = evaluate_D̄(cache, r)
    traceless_derivative = zero(cache.traceless_derivative)

    n = nstates(cache)
    for i in eachindex(derivative)
        traceless_derivative[i] = derivative[i] - D̄[i].*I(n)
    end
    return traceless_derivative
end

function evaluate_traceless_adiabatic_derivative(cache::Union{RingPolymer_QuantumModel_Cache, RingPolymer_QuantumFrictionModel_Cache}, r::AbstractArray{T,3}) where {T}
    adiabatic_derivative = evaluate_adiabatic_derivative(cache, r)
    D̄ = evaluate_D̄(cache, r)
    traceless_adiabatic_derivative = zero(cache.traceless_adiabatic_derivative)

    n = nstates(cache)
    for i in eachindex(D̄)
        cache.traceless_adiabatic_derivative[i] = adiabatic_derivative[i] - D̄[i].*I(n)
    end
    return traceless_adiabatic_derivative
end
