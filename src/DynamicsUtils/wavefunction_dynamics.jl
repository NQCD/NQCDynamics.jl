using LinearAlgebra: LinearAlgebra, mul!, lmul!, diagind, Hermitian
using LinearAlgebra.LAPACK
using NQCModels: NQCModels

function set_single_electron_derivative!(dc, c, V, v, d, tmp)
    @. dc = -im*V * c
    for I in eachindex(v)
        mul!(tmp, d[I], c)
        lmul!(v[I], tmp)
        @. dc -= tmp
    end
    return nothing
end

function propagate_wavefunction!(σfinal, σ, v, r, sim, dt)
    propagator = get_quantum_propagator(sim, v, r, dt)

    tmp1 = sim.method.tmp_matrix_complex_rect1
    tmp2 = sim.method.tmp_matrix_complex_rect2
    copy!(tmp1, σ)
    mul!(tmp2, propagator, tmp1)
    copy!(σfinal, tmp2)
end

function get_quantum_propagator(sim, v, r, dt)
    prop = sim.method.quantum_propagator
    tmp1 = sim.method.tmp_matrix_complex_square1
    tmp2 = sim.method.tmp_matrix_complex_square2

    v = DynamicsUtils.get_hopping_velocity(sim, v)
    eigenvalues = DynamicsUtils.get_hopping_eigenvalues(sim, r)
    fill!(prop, zero(eltype(prop)))
    @inbounds for i in eachindex(eigenvalues)
        prop[i,i] = eigenvalues[i]
    end

    d = DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r)
    @inbounds for i in NQCModels.mobileatoms(sim)
        for j in NQCModels.dofs(sim)
            @. prop -= 1im * d[j,i] * v[j,i]
        end
    end

    tmp1 .= Hermitian(prop)
    vals, vecs = LAPACK.syevr!('V', 'A', 'U', tmp1, 0.0, 0.0, 0, 0, 1e-5)
    
    fill!(prop, zero(eltype(prop)))
    @inbounds for i in eachindex(vals)
        prop[i,i] = exp(-1im * vals[i] * dt)
    end


    mul!(tmp2, prop, vecs')
    mul!(prop, vecs, tmp2)


    return prop
end

function get_quantum_propagator!(prop, u, p, t)

    sim, nuclei = p
    r, v = nuclei.r, nuclei.v
    v = DynamicsUtils.get_hopping_velocity(sim, v)
    eigenvalues = DynamicsUtils.get_hopping_eigenvalues(sim, r)
    fill!(prop, zero(eltype(prop)))
    @inbounds for i in eachindex(eigenvalues)
        prop[i,i] = -im * eigenvalues[i]
    end

    d = DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r)
    
    @inbounds for i in NQCModels.mobileatoms(sim)
        for j in NQCModels.dofs(sim)
            @. prop -= d[j,i] * v[j,i]
        end
    end
    return nothing
end


