#Definitions of update functions

function update_potential!(cache::Abstract_ClassicalModel_Cache, r::AbstractMatrix)
    cache.potential .= hcat(NQCModels.potential(cache.model, r))
end

function update_potential!(cache::Abstract_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    @views @inbounds for i in beads(cache)
        cache.potential[i] .= hcat(NQCModels.potential(cache.model, r[:,:,i]))
    end
    return nothing
end

function update_potential!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    NQCModels.potential!(cache.model, cache.potential, r)
end

function update_potential!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    @views @inbounds for i in beads(cache)
        NQCModels.potential!(cache.model, cache.potential[i], r[:,:,i])
    end
    return nothing
end

function update_derivative!(cache::Abstract_ClassicalModel_Cache, r::AbstractMatrix)
    NQCModels.derivative!(cache.model, cache.derivative, r::AbstractMatrix)
end

function update_derivative!(cache::Abstract_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    @views @inbounds for i in beads(cache)
        NQCModels.derivative!(cache.model, cache.derivative[:,:,i], r[:,:,i])
    end
    return nothing
end

function update_derivative!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    NQCModels.derivative!(cache.model, cache.derivative, r::AbstractMatrix)
end

function update_derivative!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    @views @inbounds for i in beads(cache)
        NQCModels.derivative!(cache.model, cache.derivative[:,:,i], r[:,:,i])
    end
    return nothing
end

function update_friction!(cache::Abstract_ClassicalModel_Cache, R::AbstractMatrix)
    NQCModels.friction!(cache.model, cache.friction, R::AbstractMatrix)
end

function update_friction!(cache::Abstract_ClassicalModel_Cache, R::AbstractArray{T,3}) where {T}
    @views for i in axes(R, 3)
        NQCModels.friction!(cache.model, cache.friction[:,:,i], R[:,:,i])
    end
    return nothing
end

function update_friction!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    μ = NQCModels.fermilevel(cache.model)

    if cache.friction_method isa WideBandExact
        potential = get_potential(cache, r)
        derivative = get_derivative(cache, r)

        fill_friction_tensor!(cache.friction, cache.friction_method, potential, derivative, r, μ)

    elseif cache.friction_method isa Nothing
        cache.friction .= zero(cache.friction)

    else
        ∂H = get_adiabatic_derivative(cache, r)
        eigen = get_eigen(cache, r)

        fill_friction_tensor!(cache.friction, cache.friction_method, eigen, ∂H, r, μ)
    end
end

function update_friction!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    μ = NQCModels.fermilevel(cache.model)
    if cache.friction_method isa WideBandExact
        potential = get_potential(cache, r)
        derivative = get_derivative(cache, r)

        @views for i in axes(r, 3)
            fill_friction_tensor!(cache.friction[:,:,i], cache.friction_method, potential[i], derivative[:,:,i], r[:,:,i], μ)
        end
    elseif cache.friction_method isa Nothing
        cache.friction .= zero(cache.friction)

    else
        ∂H = get_adiabatic_derivative(cache, r)
        eigen = get_eigen(cache, r)

        @views for i in axes(r, 3)
            fill_friction_tensor!(cache.friction[:,:,i], cache.friction_method, eigen[i], ∂H[:,:,i], r[:,:,i], μ)
        end
    end
end

function update_centroid_friction!(cache::Abstract_ClassicalModel_Cache, R::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(R)
    NQCModels.friction!(cache.model, cache.centroid_friction, centroid)
end

function update_centroid_friction!(cache::Abstract_QuantumModel_Cache, R::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(R)
    μ = NQCModels.fermilevel(cache.model)

    if cache.friction_method isa WideBandExact
        potential = get_centroid_potential(cache, centroid)
        derivative = get_centroid_derivative(cache, centroid)

        fill_friction_tensor!(cache.centroid_friction, cache.friction_method, potential, derivative, centroid, μ)

    elseif cache.friction_method isa Nothing
        cache.centroid_friction .= zero(cache.centroid_friction)
    else
        ∂H = get_centroid_adiabatic_derivative(cache, R)
        eigen = get_centroid_eigen(cache, R)

        fill_friction_tensor!(cache.centroid_friction, cache.friction_method, ∂H, eigen, centroid, μ)
    end
end

function update_eigen!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)

    # Ensure symmetry explicitly
    cache.tmp_mat .= cache.potential

    FastLapackInterface.syevr!(cache.eigen, 'V', 'A', 'U', cache.tmp_mat, 0.0, 0.0, 0, 0, 1e-12)
    correct_phase!(cache, cache.eigen)
    return nothing
end

function update_eigen!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    potential = get_potential(cache, r)

    @inbounds for i in beads(cache)
        cache.tmp_mat .= potential[i]
        FastLapackInterface.syevr!(cache.eigen[i], 'V', 'A', 'U', cache.tmp_mat, 0.0, 0.0, 0, 0, 1e-12)
        correct_phase!(cache.phase_ref[i], cache.eigen[i])
    end
    return nothing
end

function update_adiabatic_derivative!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    U = get_eigen(cache, r).Z
    diabatic_derivative = get_derivative(cache, r)

    for I in eachindex(diabatic_derivative)
        cache.adiabatic_derivative[I] .= U' * diabatic_derivative[I] * U
    end
    return nothing
end

function update_adiabatic_derivative!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    derivative = get_derivative(cache, r)
    eigen = get_eigen(cache, r)

    for i in axes(derivative, 3) # Beads
        for j in axes(derivative, 2) # Atoms
            for k in axes(derivative, 1) # DoFs
                cache.adiabatic_derivative[k,j,i] .= eigen[i].Z' * derivative[k,j,i] * eigen[i].Z
            end
        end
    end
    return nothing
end

function update_centroid_adiabatic_derivative!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T, 3}) where {T}
    centroid_derivative = get_centroid_derivative(cache, r)
    centroid_eigen = get_centroid_eigen(cache, r)
    for I in eachindex(centroid_derivative)
        cache.centroid_adiabatic_derivative[I] .= centroid_eigen.Z' * centroid_derivative[I] * centroid_eigen.Z
    end
    return nothing
end

function update_inverse_difference_matrix!(out, eigenvalues)
    @inbounds for i in eachindex(eigenvalues)
        for j in eachindex(eigenvalues)
            # Checking if the eigenvalues are different, otherwise set them to 0. Blame Henry for this!
            out[j,i] = 1 / (eigenvalues[i] - eigenvalues[j]) .* (eigenvalues[i] .!= eigenvalues[j]) 
        end
    end
    return nothing
end

"""
# References

- HammesSchifferTully_JChemPhys_101_4657_1994 Eq. (32)
- SubotnikBellonzi_AnnuRevPhyschem_67_387_2016, section 2.3
"""
function update_nonadiabatic_coupling!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    eigen = get_eigen(cache, r)
    adiabatic_derivative = get_adiabatic_derivative(cache, r)

    update_inverse_difference_matrix!(cache.tmp_mat, eigen.w)

    nonadiabatic_coupling_loop!(cache, adiabatic_derivative, cache.model)
    
    return nothing
end

function nonadiabatic_coupling_loop!(cache, adiabatic_derivative, model)
    @inbounds for I in NQCModels.dofs(cache)
        @. cache.nonadiabatic_coupling[I] = adiabatic_derivative[I] * cache.tmp_mat
    end
end

function nonadiabatic_coupling_loop!(cache, adiabatic_derivative, model::NQCModels.SpinBoson)
    @inbounds for N in eachindex(model.cⱼ)
        @inbounds for I in NQCModels.dofs(cache)
            @. cache.nonadiabatic_coupling[N] = adiabatic_derivative[N] * cache.tmp_mat
        end
    end
end

function update_nonadiabatic_coupling!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    eigen = get_eigen(cache, r)
    adiabatic_derivative = get_adiabatic_derivative(cache, r)

    @inbounds for i in beads(cache)
        update_inverse_difference_matrix!(cache.tmp_mat, eigen[i].w)
        for j in mobileatoms(cache)
            for k in dofs(cache)
                @. cache.nonadiabatic_coupling[k,j,i] = adiabatic_derivative[k,j,i] * cache.tmp_mat
            end
        end
    end
    return nothing
end

function update_centroid_nonadiabatic_coupling!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    adiabatic_derivative = get_centroid_adiabatic_derivative(cache, r)
    eigen = get_centroid_eigen(cache, r)
    update_inverse_difference_matrix!(cache.tmp_mat, eigen.w)

    @inbounds for j in mobileatoms(cache)
        for k in dofs(cache)
            @. cache.centroid_nonadiabatic_coupling[j,k] = adiabatic_derivative[j,k] * cache.tmp_mat
        end
    end
    return nothing
end

#RingPolymer specific functions

#= function update_centroid!(cache::Abstract_Cache, r::AbstractArray{T,3}) where {T}
    cache.centroid .= RingPolymerArrays.get_centroid(r)
end =#

function update_centroid_potential!(cache::Abstract_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(r)
    cache.centroid_potential .= hcat(NQCModels.potential(cache.model, centroid))
    return nothing
end

function update_centroid_potential!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(r)
    NQCModels.potential!(cache.model, cache.centroid_potential, centroid)
    return nothing
end

function update_centroid_derivative!(cache::Abstract_Cache, r::AbstractArray{T,3}) where {T}
    centroid = RingPolymerArrays.get_centroid(r)
    NQCModels.derivative!(cache.model, cache.centroid_derivative, centroid)
    return nothing
end

function update_centroid_eigen!(cache::Abstract_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    potential = get_centroid_potential(cache, r)
    cache.tmp_mat .= potential
    FastLapackInterface.syevr!(cache.centroid_eigen, 'V', 'A', 'U', cache.tmp_mat, 0.0, 0.0, 0, 0, 1e-12)
    correct_phase!(cache.centroid_phase_ref, cache.centroid_eigen)
    return nothing
end

#update tilde quantities

function update_V̄!(cache::RingPolymer_QuantumModel_Cache, r)
    potential = get_potential(cache, r)
    for i in 1:length(cache.V̄)
        cache.V̄[i] = tr(potential[i]) / nstates(cache.model)
    end
    return nothing
end

function update_V̄!(cache::RingPolymer_QuantumFrictionModel_Cache, r)
    potential = get_potential(cache, r)
    for i in 1:length(cache.V̄)
        cache.V̄[i] = tr(potential[i]) / nstates(cache.model)
    end
    return nothing
end

function update_D̄!(cache::RingPolymer_QuantumModel_Cache, r)
    derivative = get_derivative(cache, r)
    for I in eachindex(derivative)
        cache.D̄[I] = tr(derivative[I]) / nstates(cache.model)
    end
    return nothing
end

function update_D̄!(cache::RingPolymer_QuantumFrictionModel_Cache, r)
    derivative = get_derivative(cache, r)
    for I in eachindex(derivative)
        cache.D̄[I] = tr(derivative[I]) / nstates(cache.model)
    end
    return nothing
end

#update traceless quantities

function update_traceless_potential!(cache::RingPolymer_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    n = nstates(cache)
    potential = get_potential(cache, r)
    V̄ = get_V̄(cache, r)
    for i in eachindex(potential)
        cache.traceless_potential[i] = potential[i] - V̄[i].*I(n)
    end
    return nothing
end

function update_traceless_potential!(cache::RingPolymer_QuantumFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    n = nstates(cache.model)
    potential = get_potential(cache, r)
    V̄ = get_V̄(cache, r)
    for i in eachindex(potential)
        cache.traceless_potential[i] = potential[i] - V̄[i].*I(n)
    end
    return nothing
end

function update_traceless_derivative!(cache::RingPolymer_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    n = nstates(cache.model)
    derivative = get_derivative(cache, r)
    D̄ = get_D̄(cache, r)
    for i in eachindex(derivative)
        cache.traceless_derivative[i] = derivative[i] - D̄[i].*I(n)
    end
    return nothing
end

function update_traceless_derivative!(cache::RingPolymer_QuantumFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    n = nstates(cache.model)
    derivative = get_derivative(cache, r)
    D̄ = get_D̄(cache, r)
    for i in eachindex(derivative)
        cache.traceless_derivative[i] = derivative[i] - D̄[i].*I(n)
    end
    return nothing
end

function update_traceless_adiabatic_derivative!(cache::RingPolymer_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    n = nstates(cache.model)
    adiabatic_derivative = get_adiabatic_derivative(cache, r)
    D̄ = get_D̄(cache, r)
    for i in eachindex(D̄)
        cache.traceless_adiabatic_derivative[i] = adiabatic_derivative[i] - D̄[i].*I(n)
    end
    return nothing
end

function update_traceless_adiabatic_derivative!(cache::RingPolymer_QuantumFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    n = nstates(cache.model)
    adiabatic_derivative = get_adiabatic_derivative(cache, r)
    D̄ = get_D̄(cache, r)
    for i in eachindex(D̄)
        cache.traceless_adiabatic_derivative[i] = adiabatic_derivative[i] - D̄[i].*I(n)
    end
    return nothing
end

"""
updates all model properties stored in the cache for the current position `r`.
# Properties that may be updated:
- Diabatic potential
- Diabatic derivative
- Eigenvalues and eigenvectors
- Adiabatic derivative
- Nonadiabatic coupling
- Friction tensor

- Centroid equivalents of the above 
"""
function update_cache!(cache::Abstract_ClassicalModel_Cache, r::AbstractMatrix)
    update_potential!(cache, r)
    update_derivative!(cache, r)
    return nothing
end

function update_cache!(cache::ClassicalFrictionModel_Cache, r::AbstractMatrix)
    update_potential!(cache, r)
    update_derivative!(cache, r)
    update_friction!(cache, r)
    return nothing
end

function update_cache!(cache::RingPolymer_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    update_potential!(cache, r)
    update_derivative!(cache, r)

    update_centroid!(cache, r)
    return nothing
end

function update_cache!(cache::RingPolymer_ClassicalFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    update_potential!(cache, r)
    update_derivative!(cache, r)
    update_friction!(cache, r)

    update_centroid!(cache, r)
    return nothing
end

function update_cache!(cache::Abstract_QuantumModel_Cache, r::AbstractMatrix)
    # println("before:")
    # println(cache)
    update_potential!(cache, r)
    update_derivative!(cache, r)
    update_eigen!(cache, r)
    update_adiabatic_derivative!(cache, r)
    update_nonadiabatic_coupling!(cache, r)
    # println("after:")
    # println(cache)
    return nothing
end

function update_cache!(cache::QuantumFrictionModel_Cache, r::AbstractMatrix)
    update_potential!(cache, r)
    update_derivative!(cache, r)
    update_eigen!(cache, r)
    update_adiabatic_derivative!(cache, r)
    update_nonadiabatic_coupling!(cache, r)
    update_friction!(cache, r)
    return nothing
end

function update_cache!(cache::RingPolymer_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    update_potential!(cache, r)
    update_derivative!(cache, r)
    update_eigen!(cache, r)
    update_adiabatic_derivative!(cache, r)
    update_nonadiabatic_coupling!(cache, r)

    update_V̄!(cache, r)
    update_D̄!(cache, r)

    update_traceless_potential!(cache,r)
    update_traceless_derivative!(cache,r)
    update_traceless_adiabatic_derivative!(cache,r)

    update_centroid!(cache, r)
    return nothing
end

function update_cache!(cache::RingPolymer_QuantumFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    update_potential!(cache, r)
    update_derivative!(cache, r)
    update_eigen!(cache, r)
    update_adiabatic_derivative!(cache, r)
    update_nonadiabatic_coupling!(cache, r)
    update_friction!(cache, r)

    update_V̄!(cache, r)
    update_D̄!(cache, r)

    update_traceless_potential!(cache,r)
    update_traceless_derivative!(cache,r)
    update_traceless_adiabatic_derivative!(cache,r)

    update_centroid!(cache, r)
    return nothing
end

"""
updates all model properties stored in the cache for the current centroid position `r_centroid`.
# Properties that may be updated:
- Diabatic potential
- Diabatic derivative
- Eigenvalues and eigenvectors
- Adiabatic derivative
- Nonadiabatic coupling
- Friction tensor 
"""
function update_centroid!(cache::RingPolymer_ClassicalModel_Cache, r::AbstractArray{T,3}) where {T}
    cache.centroid .= RingPolymerArrays.get_centroid(r)
    update_centroid_potential!(cache, r)
    update_centroid_derivative!(cache, r)
    return nothing
end

function update_centroid!(cache::RingPolymer_ClassicalFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    cache.centroid .= RingPolymerArrays.get_centroid(r)
    update_centroid_potential!(cache, r)
    update_centroid_derivative!(cache, r)
    update_centroid_friction!(cache, r)
    return nothing
end

function update_centroid!(cache::RingPolymer_QuantumModel_Cache, r::AbstractArray{T,3}) where {T}
    cache.centroid .= RingPolymerArrays.get_centroid(r)
    update_centroid_potential!(cache, r)
    update_centroid_derivative!(cache, r)
    update_centroid_eigen!(cache, r)
    update_centroid_adiabatic_derivative!(cache, r)
    update_centroid_nonadiabatic_coupling!(cache, r)
    return nothing
end

function update_centroid!(cache::RingPolymer_QuantumFrictionModel_Cache, r::AbstractArray{T,3}) where {T}
    cache.centroid .= RingPolymerArrays.get_centroid(r)
    update_centroid_potential!(cache, r)
    update_centroid_derivative!(cache, r)
    update_centroid_eigen!(cache, r)
    update_centroid_adiabatic_derivative!(cache, r)
    update_centroid_nonadiabatic_coupling!(cache, r)
    update_centroid_friction!(cache, r)
    return nothing
end
