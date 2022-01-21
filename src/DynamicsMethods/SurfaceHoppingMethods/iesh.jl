using LinearAlgebra: mul!, det

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
    n_electrons::Int
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    rescaling::Symbol
    function AdiabaticIESH{T}(states::Integer, n_electrons::Integer, rescaling::Symbol) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)

        new{T}(hopping_probability, state, new_state, proposed_state,n_electrons, overlap, tmp, rescaling)
    end
end

<<<<<<< HEAD
function NonadiabaticMolecularDynamics.Simulation{IESH_type}(atoms::Atoms{S,T}, model::Model; n_electrons, rescaling=:standard, kwargs...) where {S,T,IESH_type<:AbstractIESH}
    NonadiabaticMolecularDynamics.Simulation(atoms, model, IESH_type{T}(NonadiabaticModels.nstates(model), n_electrons, rescaling); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:AdiabaticIESH}, v, r)
    ψ = zeros(NonadiabaticModels.nstates(sim.calculator.model), sim.method.n_electrons)
    
=======
function NQCDynamics.Simulation{IESH}(atoms::Atoms{S,T}, model::Model; n_electrons, kwargs...) where {S,T}
    NQCDynamics.Simulation(atoms, model, IESH{T}(NQCModels.nstates(model), n_electrons); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:IESH}, v, r)
    ψ = zeros(NQCModels.nstates(sim.calculator.model), sim.method.n_electrons)

>>>>>>> master
    for i=1:sim.method.n_electrons
        ψ[i,i] = 1
    end
    state = collect(1:sim.method.n_electrons)

    SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=ψ, σimag=zero(ψ)), state)
end

"""
Set the acceleration due to the force from the currently occupied states.
See Eq. 12 of Shenvi, Tully JCP 2009 paper.
"""
function acceleration!(dv, v, r, sim::Simulation{<:AbstractIESH}, t, state)
    dv .= zero(eltype(dv))
    NonadiabaticModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)
    for I in eachindex(dv)
        for k in state
            # Contribution to the force from each occupied state `k`
            dv[I] -= sim.calculator.adiabatic_derivative[I][k, k]
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
    @views for i in axes(dσ, 2)       
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
    hamiltonian = LinearAlgebra.diagm(sim.calculator.eigen.values)
    hamiltonian = complex(hamiltonian)
    for I in eachindex(v)
        hamiltonian .-= 1im .* sim.calculator.nonadiabatic_coupling[I] .* v[I]
    end
    vals, vecs = LinearAlgebra.eigen(hamiltonian)
    propagator = LinearAlgebra.diagm(exp.(-im .* vals .* dt))
    propagator = vecs * propagator * vecs'
    return propagator
end

"""
Hopping probability according to equation 21 in Shenvi, Roy, Tully 2009.
"""
function evaluate_hopping_probability!(sim::Simulation{<:AbstractIESH}, u, dt)

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    v = DynamicsUtils.get_velocities(u)

    S = sim.method.overlap
    proposed_state = sim.method.proposed_state
    prob = sim.method.hopping_probability
    d = sim.calculator.nonadiabatic_coupling

    compute_overlap!(sim, S, ψ, sim.method.state)

    det_current = det(S)
    Akk = abs2(det_current)

    fill!(prob, zero(eltype(prob)))
    for i in axes(prob, 2) # each electron
        for j in axes(prob, 1) # each basis state
            if j ∉ sim.method.state
                copyto!(proposed_state, sim.method.state)
                proposed_state[i] = j
                Akj = calculate_Akj(sim, S, ψ, det_current, proposed_state)
                for I in eachindex(v)
                    prob[j,i] += 2 * v[I] * real(Akj/Akk) * d[I][sim.method.state[i], j] * dt
                end
            end
        end
    end

    clamp!(prob, 0, 1) # Restrict probabilities between 0 and 1

    return nothing
end

"Equation 17 in Shenvi, Roy, Tully 2009. Uses equations 19 and 20."
function calculate_Akj(sim::Simulation{<:AbstractIESH}, S::AbstractMatrix, ψ::AbstractMatrix, detS::Number, new_state::Vector)

    compute_overlap!(sim, S, ψ, new_state)
    det_new = det(S)

    return detS*conj(det_new)
end

"Equation 20 in Shenvi, Roy, Tully 2009."
function compute_overlap!(::Simulation{<:AdiabaticIESH}, S::Matrix, ψ, state)
    for i in axes(S, 2) # each electron
        for j in axes(S, 1) # each electron
            S[j,i] = ψ[state[j],i]
        end
    end
    return nothing
end

function select_new_state(sim::AbstractSimulation{<:AbstractIESH}, u)::Vector{Int}

    prob = sim.method.hopping_probability

    random = rand()
    cumulative = 0
    for i in axes(prob, 2) # each electron
        for j in axes(prob, 1) # each basis state
            if j ∉ sim.method.state
                cumulative += prob[j,i]
                if random < cumulative
                    copyto!(sim.method.proposed_state, sim.method.state)
                    sim.method.proposed_state[i] = j
                    return sim.method.proposed_state
                end
            end
        end
    end

    return sim.method.state
end

function Estimators.diabatic_population(sim::Simulation{<:AdiabaticIESH}, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigen.vectors

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    diabatic_ψ = zero(ψ)
    @views for i in axes(ψ, 2) # each electron
        diabatic_ψ[:,i] .= U*ψ[:,i]
    end
    diabatic_population = sum(abs2, diabatic_ψ, dims=2)

    return diabatic_population
end

<<<<<<< HEAD
function Estimators.adiabatic_population(sim::Simulation{<:AdiabaticIESH}, u)
    population = zeros(NonadiabaticModels.nstates(sim.calculator.model))
    for i=1:sim.method.n_electrons
        population .+= abs2.(DynamicsUtils.get_quantum_subsystem(u)[:,i])
    end
=======
function Estimators.adiabatic_population(sim::Simulation{<:IESH}, u)
    population = zeros(NQCModels.nstates(sim.calculator.model))
    population[u.state] .= 1
>>>>>>> master
    return population
end

unpack_states(sim::Simulation{<:AbstractIESH}) = symdiff(sim.method.new_state, sim.method.state)

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:AbstractIESH}, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    potential = NonadiabaticModels.state_independent_potential(sim.calculator.model, DynamicsUtils.get_positions(u))
    for i in u.state
        potential += sim.calculator.eigen.values[i]
    end
    return potential
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:AbstractIESH}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    potential = DynamicsUtils.classical_potential_energy(sim, u)
    return kinetic + potential
end
