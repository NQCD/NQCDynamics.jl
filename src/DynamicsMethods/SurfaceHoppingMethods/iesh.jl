using LinearAlgebra: mul!, det

export IESH

"""
    IESH{T} <: SurfaceHopping

Independent electron surface hopping.

# References
- [Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009)](https://doi.org/10.1063/1.3125436)
- [Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009)](https://doi.org/10.1063/1.3122989)

"""
struct IESH{T} <: SurfaceHopping
    hopping_probability::Matrix{T}
    state::Vector{Int}
    new_state::Vector{Int}
    proposed_state::Vector{Int}
    n_electrons::Int
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    function IESH{T}(states::Integer, n_electrons::Integer) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)

        new{T}(hopping_probability, state, new_state, proposed_state,n_electrons, overlap, tmp)
    end
end

function NonadiabaticMolecularDynamics.Simulation{IESH}(atoms::Atoms{S,T}, model::Model; n_electrons, kwargs...) where {S,T}
    NonadiabaticMolecularDynamics.Simulation(atoms, model, IESH{T}(NonadiabaticModels.nstates(model), n_electrons); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:IESH}, v, r)
    ψ = zeros(NonadiabaticModels.nstates(sim.calculator.model), sim.method.n_electrons)

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
function acceleration!(dv, v, r, sim::Simulation{<:IESH}, t, state)
    dv .= zero(eltype(dv))
    for I in eachindex(dv)
        for k in state
            # Contribution to the force from each occupied state `k`
            # Goes to Calculator.jl
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
function set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:IESH})
    V = sim.calculator.eigenvalues
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

"""
Hopping probability according to equation 21 in Shenvi, Roy, Tully 2009.
"""
function evaluate_hopping_probability!(sim::Simulation{<:IESH}, u, dt)

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    v = DynamicsUtils.get_velocities(u)

    S = sim.method.overlap
    proposed_state = sim.method.proposed_state
    prob = sim.method.hopping_probability
    d = sim.calculator.nonadiabatic_coupling

    compute_overlap!(S, ψ, sim.method.state)
    det_current = det(S)
    Akk = abs2(det_current)

    fill!(prob, zero(eltype(prob)))
    for i in axes(prob, 2) # each electron
        for j in axes(prob, 1) # each basis state
            if j ∉ sim.method.state
                copyto!(proposed_state, sim.method.state)
                proposed_state[i] = j
                Akj = calculate_Akj(S, ψ, det_current, proposed_state)
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
function calculate_Akj(S, ψ, detS, new_state)

    compute_overlap!(S, ψ, new_state)
    det_new = det(S)

    return detS*conj(det_new)
end

"Equation 20 in Shenvi, Roy, Tully 2009."
function compute_overlap!(S::Matrix, ψ, state)
    for i in axes(S, 2) # each electron
        for j in axes(S, 1) # each electron
            S[j,i] = ψ[state[j],i]
        end
    end
    return nothing
end

function select_new_state(sim::AbstractSimulation{<:IESH}, u)::Vector{Int}

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

function Estimators.diabatic_population(sim::Simulation{<:IESH}, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigenvectors

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    diabatic_ψ = zero(ψ)
    @views for i in axes(ψ, 2) # each electron
        diabatic_ψ[:,i] .= U*ψ[:,i]
    end
    diabatic_population = sum(abs2.(diabatic_ψ), dims=2)

    return diabatic_population
end

function Estimators.adiabatic_population(sim::Simulation{<:IESH}, u)
    population = zeros(NonadiabaticModels.nstates(sim.calculator.model))
    population[u.state] .= 1
    return population
end

# function rescale_velocity!(sim::Simulation{<:IESH}, u)::Bool
#     # ShakibHuo_JPhysChemLett_8_3073_2017_rescale
#     new_state, old_state = symdiff(sim.method.new_state, sim.method.state)
#     velocity = DynamicsUtils.get_velocities(u)
    
#     # Calculate difference in eigenvalues between old and new state, weighed by mass:
#     # c = (E_{new} - E_{old})/mass
#     c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
#     # a = d^2_{new,old}/mass (where d is the nonadiabatic coupling)
#     # b = v*d_{new,old}
#     a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
#     discriminant = b.^2 .- 2a.*c

#     # If smaller zero, hop is frustrated
#     any(discriminant .< 0) && return false

#     root = sqrt.(discriminant)
#     plus = (b .+ root) ./ a
#     minus = (b .- root) ./ a 
#     velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus
#     perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

#     return true
# end

# function calculate_potential_energy_change(calc::AbstractDiabaticCalculator, new_state::Integer, current_state::Integer)
#     return calc.eigenvalues[new_state] - calc.eigenvalues[current_state]
# end

# function evaluate_a_and_b(sim::Simulation{<:SurfaceHopping}, velocity, new_state, old_state)
#     a = zeros(length(sim.atoms))
#     b = zero(a)
#     @views for i in range(sim.atoms)
#         coupling = [sim.calculator.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
#         a[i] = dot(coupling, coupling) / sim.atoms.masses[i]
#         b[i] = dot(velocity[:,i], coupling)
#     end
#     return (a, b)
# end

# function perform_rescaling!(sim::Simulation{<:SurfaceHopping}, velocity, velocity_rescale, new_state, old_state)
#     for i in range(sim.atoms)
#         coupling = [sim.calculator.nonadiabatic_coupling[j,i][new_state, old_state] for j=1:sim.DoFs]
#         velocity[:,i] .-= velocity_rescale[i] .* coupling ./ sim.atoms.masses[i]
#     end
#     return nothing
# end