
export DiabaticIESH

struct DiabaticIESH{T} <: AbstractIESH
    hopping_probability::Matrix{T}
    state::Vector{Int}
    new_state::Vector{Int}
    proposed_state::Vector{Int}
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    tmp_matrix_complex1::Matrix{Complex{T}}
    tmp_matrix_complex2::Matrix{Complex{T}}
    rescaling::Symbol
    tmp_vector_int::Vector{Int}
    v_dot_d::Matrix{T}
    unoccupied::Vector{Int}
    function DiabaticIESH{T}(states::Integer, n_electrons::Integer, rescaling::Symbol) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        tmp_matrix_complex1 = zeros(Complex{T}, states, n_electrons)
        tmp_matrix_complex2 = zeros(Complex{T}, states, n_electrons)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)
        tmp_vector_int = zeros(Int, n_electrons)
        v_dot_d = zeros(T, states, n_electrons)
        unoccupied = zeros(Int, states - n_electrons)

        new{T}(hopping_probability, state, new_state, proposed_state, overlap,
            tmp, tmp_matrix_complex1, tmp_matrix_complex2, rescaling, tmp_vector_int, v_dot_d, unoccupied)
    end
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:DiabaticIESH}, v, r)

    ψ = zeros(NQCModels.nstates(sim), NQCModels.nelectrons(sim))
    for i in eachelectron(sim)
        ψ[i,i] = 1
    end
    state = collect(eachelectron(sim))

    Calculators.update_electronics!(sim.calculator, r)
    U = sim.calculator.eigen.vectors
    @views for i in eachelectron(sim)
        ψ[:,i] .= U * ψ[:,i]
    end
    # ψ = U[:,state] The same as this, transforming population from adiabatic to diabatic

    SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=ψ, σimag=zero(ψ)), state)
end

function propagate_wavefunction!(σfinal, σ, v, sim::Simulation{<:DiabaticIESH}, dt)
    U = sim.calculator.eigen.vectors
    converted = U'*σ
    propagator = exp.(-im .* sim.calculator.eigen.values .* dt)
    for i in axes(σfinal, 2)
        for j in axes(σfinal, 1)
            σfinal[j,i] = propagator[j] * converted[j,i]
        end
    end
    σfinal .= U * σfinal
end

function DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:DiabaticIESH})
    V = sim.calculator.potential
    for (dc, c) in zip(eachcol(dσ), eachcol(σ))
        mul!(dc, V, c)
        lmul!(-im, dc)
    end
end

"Equation 20 in Shenvi, Roy, Tully 2009."
function compute_overlap!(sim::Simulation{<:DiabaticIESH}, S::Matrix, ψ, state)
    tmp1 = sim.method.tmp_matrix_complex1
    tmp2 = sim.method.tmp_matrix_complex2
    ϕ = @view sim.calculator.eigen.vectors[:,state]
    copy!(tmp1, ϕ) # Copy to complex array
    copy!(tmp2, ψ) # Copy to complex array
    LinearAlgebra.mul!(S, tmp1', tmp2) # Using tmp arrays for performance
end

function Estimators.adiabatic_population(sim::Simulation{<:DiabaticIESH}, u)
    population = zeros(NQCModels.nstates(sim.calculator.model))
    Calculators.update_electronics!(sim.calculator, DynamicsUtils.get_positions(u))
    U = sim.calculator.eigen.vectors
    for i in eachelectron(sim)
        adiabatic = U' * DynamicsUtils.get_quantum_subsystem(u)[:,i]
        population .+= abs2.(adiabatic)
    end
    return population
end

function Estimators.diabatic_population(::Simulation{<:DiabaticIESH}, u)
    return sum(abs2, DynamicsUtils.get_quantum_subsystem(u); dims=2)
end
