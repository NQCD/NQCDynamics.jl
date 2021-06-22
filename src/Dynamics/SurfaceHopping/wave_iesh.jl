
export IESH
export SurfaceHoppingVariablesIESH

"""This module controles how IESH is executed. For a description of IESH, see e.g.
Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009) and 
Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009).

The density matrix is set up in surface_hopping_variables.jl
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

function SurfaceHoppingVariablesIESH(v::AbstractArray, r::AbstractArray, n_states::Integer, n_electrons::Integer)
    ψ = zeros(Complex{eltype(r)}, n_states, n_electrons)

    for i=1:n_electrons
        ψ[i,i] = 1
    end
    state = collect(1:n_electrons)

    SurfaceHoppingVariables(ArrayPartition(v, r, ψ), state)
end

"""This is part of the adiabatic propagation of the nuclei.
   See Eq. 12 of Shenvi, Tully JCP 2009 paper."""
function acceleration!(dv, v, r, sim::Simulation{<:IESH}, t, state)
    dv .= 0.0
    for i in axes(dv, 2)
        for j in axes(dv, 1)        
            for k in state
                # Contribution to the force from each occupied state `k`
                dv[j,i] -= sim.calculator.adiabatic_derivative[j,i][k, k] / sim.atoms.masses[i]
            end
        end
    end
    return nothing
end

"""Propagation of electronic wave function happens according to Eq. (14) 
   in the Shenvi, Tully paper (JCP 2009)
   The extended formula is taken from HammesSchifferTully_JChemPhys_101_4657_1994, Eq. (16):
   iħ d ψ_{j}/dt = ∑_j ψ_{m}(V_{jm} - i v d_{jm})
   Where v is the velocity. See also Tully_JChemPhys_93_1061_1990, Eq. (7). 
   For the IESH case, since each electron is treated independently, the equation
   above needs to be evaluated for each electron, so it becomes:
   iħ d ψ^{K}_j/dt =  V_{j,j}ψ^{K}_j - i v ∑_m d_{m,j}*ψ^{K}_m)
   K is the number of the electron and 
   j, m run over the number of states
   """
function set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:IESH})
    V = sim.calculator.eigenvalues
    d = sim.calculator.nonadiabatic_coupling
    @views for i in axes(dσ, 2)       
        set_single_electron_derivative!(dσ[:,i], σ[:,i], V, v, d, sim.method.tmp)
    end
end

function set_single_electron_derivative!(dc, c, V, v, d, tmp)
    @. dc = -im*V * c
    for I in eachindex(v)
        mul!(tmp, d[I], c)
        lmul!(v[I], tmp)
        @. dc -= tmp
    end
    return nothing
end

"""Hopping probability according to Eq.s (21), (17), (18) in Shenvi/Roy/Tully JCP 2009 paper.
   The density matrix is used there:
   σ_{i,j} = ∑_K c_{K,i}*c_{K,j}^*
   K is the index of the independent electrons
   i,j are the electronic states.
   The equation for the hopping probability is:
   g_{k,j} = Max(\frac{-2 Real(σ_{k,j} v d_{j,k})}{σ_{k,k}})
   """
function evaluate_hopping_probability!(sim::Simulation{<:IESH}, u, dt)

    ψ = get_quantum_subsystem(u)
    v = get_velocities(u)

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

function get_diabatic_population(sim::Simulation{<:IESH}, u)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigenvectors

    ψ = get_quantum_subsystem(u)
    diabatic_ψ = zero(ψ)
    @views for i in axes(ψ, 2) # each electron
        diabatic_ψ[:,i] .= U*ψ[:,i]
    end
    diabatic_population = sum(abs2.(diabatic_ψ), dims=2)

    return diabatic_population
end

function get_adiabatic_population(sim::Simulation{<:IESH}, u)
    population = zeros(sim.calculator.model.n_states)
    population[u.state] .= 1
    return population
end

function rescale_velocity!(sim::Simulation{<:IESH}, u)::Bool
    new_state, old_state = symdiff(sim.method.new_state, sim.method.state)
    velocity = get_velocities(u)
    
    c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
    a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
    discriminant = b.^2 .- 2a.*c

    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    plus = (b .+ root) ./ a
    minus = (b .- root) ./ a 
    velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus
    perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

    return true
end

# function impurity_summary(model::DiabaticModel, R::AbstractMatrix, state::AbstractArray, σ::AbstractArray)
#     """Calculate impurity population according to MiaoSubotnik_JChemPhys_150_041711_2019"""

#     eig_vec = zeros(model.n_states,model.n_states)

#     eival = zeros(model.n_states)
#     σad = zeros(Complex, model.n_states, model.n_states)
#     tmp = 0
#     σdia = zeros(Complex, model.n_states, model.n_states)
#     eig_array = zeros(4+model.n_states)
#     V = Hermitian(zeros(model.n_states,model.n_states))
#     dvect = zeros(model.n_states)
#     dvect[2] = 1
    
#     # get potential
#     potential!(model,V,R)
#     eig_vec .= eigvecs(V)
#     ieig = inv(eig_vec)
#     # Get density matrix matrix
#     for i = 1:length(state)
#         eig_array[4 + i] = norm(σ[i,:])
#         for j = 1:length(state)
#             for k in axes(σ,2)
#                 σad[i,j] = σad[i,j] + σ[i,k]*σ[j,k]
#             end
#         end
#     end


#     #Turn into diabatic matrix for impurity population
#     σdia .= eig_vec *σad * ieig
#     eig_array[4] = real(σdia[2,2])^2 + imag(σdia[2,2])^2

#     # eig_array[4] = 0
#     # for k in axes(σ,2)
#     #     for i = 1:length(state)
#     #         σad[i,i] = σ[i,k]
#     #     end
#     #     σdia .= 0
#     #     σdia .= eig_vec *σad * ieig
#     #     eig_array[4] = eig_array[4] + real(σdia[2,2])^2 + imag(σdia[2,2])^2
#     # end
    

#     # Get the eigenvectors and values
#     eival .= eigvals(V)

#     # save position
#     eig_array[1] = R[1]
#     for i = 1:length(state)
#         # Energy
#         eig_array[2] = eig_array[2] + eival[state[i]]
#         # Hopping prob. by hopping array
#         eig_array[3] = eig_array[3] + state[i]
        
#     end

#     # # over electrons. This is adiabatic.
#     # for i in axes(σ, 2)
#     #     eig_array[4] = eig_array[4] + real(σ[2,i])^2 + imag(σ[2,i])^2
#     # end

#     # Export an array of eigenvalues with last two elements being hopping prob
#     eig_array = eig_array
# end
