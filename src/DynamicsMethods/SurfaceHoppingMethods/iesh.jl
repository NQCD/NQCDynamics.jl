using LinearAlgebra: mul!, det, Hermitian
using Unitful, UnitfulAtomic
using MuladdMacro: @muladd
using NQCDynamics: FastDeterminant
using NQCModels: eachelectron, eachstate, mobileatoms, dofs

export AdiabaticIESH

abstract type AbstractIESH <: SurfaceHopping end

"""
    IESH{T} <: SurfaceHopping

Independent electron surface hopping.

# References
- [Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009)](https://doi.org/10.1063/1.3125436)
- [Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009)](https://doi.org/10.1063/1.3122989)

"""
struct AdiabaticIESH{T} <: AbstractIESH
    hopping_probability::Matrix{T}
    state::Vector{Int}
    new_state::Vector{Int}
    proposed_state::Vector{Int}
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    rescaling::Symbol
    quantum_propagator::Matrix{Complex{T}}
    tmp_matrix_complex_square1::Matrix{Complex{T}}
    tmp_matrix_complex_square2::Matrix{Complex{T}}
    tmp_vector_int::Vector{Int}
    v_dot_d::Matrix{T}
    unoccupied::Vector{Int}
    function AdiabaticIESH{T}(states::Integer, n_electrons::Integer, rescaling::Symbol) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)
        quantum_propagator = zeros(Complex{T}, states, states)
        tmp_matrix_complex_square1 = zeros(Complex{T}, states, states)
        tmp_matrix_complex_square2 = zeros(Complex{T}, states, states)
        tmp_vector_int = zeros(Int, n_electrons)
        v_dot_d = zeros(T, states, n_electrons)
        unoccupied = zeros(Int, states - n_electrons)

        new{T}(hopping_probability, state, new_state, proposed_state, overlap, tmp, rescaling, quantum_propagator,
            tmp_matrix_complex_square1, tmp_matrix_complex_square2, tmp_vector_int, v_dot_d, unoccupied
        )
    end
end

unoccupied_states(sim::Simulation{<:AbstractIESH}) = sim.method.unoccupied

function NQCDynamics.Simulation{IESH_type}(atoms::Atoms{T}, model::Model; rescaling=:standard, kwargs...) where {T,IESH_type<:AbstractIESH}
    NQCDynamics.Simulation(atoms, model, IESH_type{T}(NQCModels.nstates(model), NQCModels.nelectrons(model), rescaling); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:AdiabaticIESH}, v, r)
    ψ = zeros(NQCModels.nstates(sim), NQCModels.nelectrons(sim))
    
    for i in eachelectron(sim)
        ψ[i,i] = 1
    end
    state = collect(eachelectron(sim))

    SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=ψ, σimag=zero(ψ)), state)
end

function DynamicsMethods.create_problem(u0, tspan, sim::AbstractSimulation{<:AbstractIESH})
    set_state!(sim.method, u0.state)
    set_unoccupied_states!(sim)
    OrdinaryDiffEq.ODEProblem(DynamicsMethods.motion!, u0, tspan, sim;
        callback=DynamicsMethods.get_callbacks(sim))
end

"""
Set the acceleration due to the force from the currently occupied states.
See Eq. 12 of Shenvi, Tully JCP 2009 paper.
"""
function DynamicsUtils.acceleration!(dv, v, r, sim::Simulation{<:AbstractIESH}, t, state)
    dv .= zero(eltype(dv))
    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)
    @inbounds for i in mobileatoms(sim)
        for j in dofs(sim)
            for k in state
                # Contribution to the force from each occupied state `k`
                dv[j,i] -= sim.calculator.adiabatic_derivative[j,i][k, k]
            end
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end

"""
Propagation of electronic wave function happens according to Eq. (14) 
in the Shenvi, Tully paper (JCP 2009)

In IESH each electron is independent so we can loop through electrons and set the
derivative one at a time, in the standard way for FSSH.
"""
function DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:AdiabaticIESH})
    V = sim.calculator.eigen.values
    d = sim.calculator.nonadiabatic_coupling
    @views for i in eachelectron(sim)
        set_single_electron_derivative!(dσ[:,i], σ[:,i], V, v, d, sim.method.tmp)
    end
end

"""
```math
\\frac{d c_k}{dt} = -iV_{kk}c_k/\\hbar - \\sum_j v d_{kj} c_j
```
"""
function set_single_electron_derivative!(dc, c, V, v, d, tmp)
    @. dc = -im*V * c
    for I in eachindex(v)
        mul!(tmp, d[I], c)
        lmul!(v[I], tmp)
        @. dc -= tmp
    end
    return nothing
end

function propagate_wavefunction!(σfinal, σ, v, sim::Simulation{<:AdiabaticIESH}, dt)
    propagator = get_quantum_propagator(sim, v, dt)
    mul!(σfinal, propagator, σ)
end

function get_quantum_propagator(sim::Simulation{<:AdiabaticIESH}, v, dt)
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

    eigs = LinearAlgebra.eigen(Hermitian(prop))
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

"""
Hopping probability according to equation 21 in Shenvi, Roy, Tully 2009.
"""
function evaluate_hopping_probability!(sim::Simulation{<:AbstractIESH}, u, dt, random)
    ψ = DynamicsUtils.get_quantum_subsystem(u)
    v = DynamicsUtils.get_velocities(u)

    S = sim.method.overlap
    proposed_state = sim.method.proposed_state
    prob = sim.method.hopping_probability
    d = sim.calculator.nonadiabatic_coupling

    compute_overlap!(sim, S, ψ, sim.method.state)

    det_current = det(S)
    Akk = abs2(det_current)
    prefactor = 2dt / real(Akk)

    fill!(prob, zero(eltype(prob)))

    evaluate_v_dot_d!(sim, v, d)

    if estimate_maximum_probability(sim, dt) < random
        return nothing # Do not evaluate probability if maximum estimate too low
    end

    @inbounds for n in eachelectron(sim)
        for m in unoccupied_states(sim)
            copy!(proposed_state, sim.method.state)
            proposed_state[n] = m
            Akj = calculate_Akj(sim, S, ψ, det_current, proposed_state)
            @muladd prob[m,n] = prob[m,n] + prefactor * real(Akj) * sim.method.v_dot_d[m,n]
        end
    end

    clamp!(prob, 0, 1) # Restrict probabilities between 0 and 1

    return nothing
end

function evaluate_v_dot_d!(sim::Simulation{<:AbstractIESH}, v, d)
    v_dot_d = sim.method.v_dot_d
    fill!(v_dot_d, zero(eltype(v_dot_d)))
    @inbounds for i in mobileatoms(sim)
        for j in dofs(sim)
            for n in eachelectron(sim)
                current_state = sim.method.state[n]
                for m in unoccupied_states(sim)
                    sim.method.v_dot_d[m,n] -= v[j,i] * d[j,i][m, current_state]
                end
            end
        end
    end
end

function estimate_maximum_probability(sim::Simulation{<:AbstractIESH}, dt)
    return sum(sim.method.v_dot_d) * 2dt
end

"Equation 17 in Shenvi, Roy, Tully 2009. Uses equations 19 and 20."
function calculate_Akj(sim::Simulation{<:AbstractIESH}, S::AbstractMatrix, ψ::AbstractMatrix, detS::Number, new_state::Vector)
    compute_overlap!(sim, S, ψ, new_state)
    det_new = FastDeterminant.det!(S, sim.method.tmp_vector_int)

    return detS*conj(det_new)
end

"Equation 20 in Shenvi, Roy, Tully 2009."
function compute_overlap!(sim::Simulation{<:AdiabaticIESH}, S::Matrix, ψ, state)
    @inbounds for i in eachelectron(sim)
        for j in eachelectron(sim)
            S[j,i] = ψ[state[j],i]
        end
    end
    return nothing
end

function select_new_state(sim::AbstractSimulation{<:AbstractIESH}, u, random)::Vector{Int}

    prob = sim.method.hopping_probability

    cumulative = 0
    for i in eachelectron(sim)
        for j in unoccupied_states(sim)
            cumulative += prob[j,i]
            if random < cumulative
                copy!(sim.method.proposed_state, sim.method.state)
                sim.method.proposed_state[i] = j
                return sim.method.proposed_state
            end
        end
    end

    return sim.method.state
end

function Estimators.diabatic_population(sim::Simulation{<:AdiabaticIESH}, u)
    eigen = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    U = eigen.vectors

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    diabatic_ψ = zero(ψ)
    @views for i in eachelectron(sim) # each electron
        diabatic_ψ[:,i] .= U*ψ[:,i]
    end
    diabatic_population = sum(abs2, diabatic_ψ, dims=2)

    return diabatic_population
end

function Estimators.adiabatic_population(sim::Simulation{<:AdiabaticIESH}, u)
    population = zeros(NQCModels.nstates(sim.calculator.model))
    for i in eachelectron(sim)
        population .+= abs2.(DynamicsUtils.get_quantum_subsystem(u)[:,i])
    end
    return population
end

unpack_states(sim::Simulation{<:AbstractIESH}) = symdiff(sim.method.new_state, sim.method.state)

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:AbstractIESH}, u)
    eigen = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = NQCModels.state_independent_potential(sim.calculator.model, DynamicsUtils.get_positions(u))
    for i in u.state
        potential += eigen.values[i]
    end
    return potential
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:AbstractIESH}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    potential = DynamicsUtils.classical_potential_energy(sim, u)
    return kinetic + potential
end

function iesh_check_hop!(u, t, integrator)::Bool
    sim = integrator.p
    random = rand()
    evaluate_hopping_probability!(sim, u, OrdinaryDiffEq.get_proposed_dt(integrator), random)
    set_new_state!(sim.method, select_new_state(sim, u, random))
    return sim.method.new_state != sim.method.state
end

function iesh_execute_hop!(integrator)
    sim = integrator.p
    if rescale_velocity!(sim, integrator.u)
        set_state!(integrator.u, sim.method.new_state)
        set_state!(sim.method, sim.method.new_state)
        set_unoccupied_states!(sim)
    end
    return nothing
end

function set_unoccupied_states!(sim::Simulation{<:AbstractIESH})
    DynamicsUtils.set_unoccupied_states!(sim.method.unoccupied, sim.method.state)
end

const IESHCallback = DiffEqBase.DiscreteCallback(iesh_check_hop!, iesh_execute_hop!;
                                                    save_positions=(false, false))

DynamicsMethods.get_callbacks(::Simulation{<:AbstractIESH}) = IESHCallback
