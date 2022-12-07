using LinearAlgebra: LinearAlgebra, mul!, lmul!, diagind, Hermitian
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

function propagate_wavefunction!(σfinal, σ, v, sim, dt)
    propagator = get_quantum_propagator(sim, v, dt)

    tmp1 = sim.method.tmp_matrix_complex_rect1
    tmp2 = sim.method.tmp_matrix_complex_rect2
    copy!(tmp1, σ)
    mul!(tmp2, propagator, tmp1)
    copy!(σfinal, tmp2)
end

function get_quantum_propagator(sim, v, dt)
    prop = sim.method.quantum_propagator
    eigenvalues = sim.calculator.eigen.values
    fill!(prop, zero(eltype(prop)))
    @inbounds for (i, I) in zip(eachindex(eigenvalues), diagind(prop))
        prop[I] = eigenvalues[i]
    end

    d = sim.calculator.nonadiabatic_coupling
    @inbounds for i in NQCModels.mobileatoms(sim)
        for j in NQCModels.dofs(sim)
            @. prop -= 1im * d[j,i] * v[j,i]
        end
    end

    eigs = LinearAlgebra.eigen!(Hermitian(prop))
    fill!(prop, zero(eltype(prop)))
    @inbounds for (i, I) in zip(eachindex(eigs.values), diagind(prop))
        prop[I] = exp(-1im * eigs.values[i] * dt)
    end

    tmp1 = sim.method.tmp_matrix_complex_square1
    tmp2 = sim.method.tmp_matrix_complex_square2

    copy!(tmp1, eigs.vectors) # Copy real->complex for faster mul!
    mul!(tmp2, prop, tmp1')
    mul!(prop, tmp1, tmp2)

    return prop
end
