using LinearAlgebra: mul!, det
using Unitful, UnitfulAtomic
using NQCModels
export IESH_Tully

    """
    IESH_Tully{T} <: SurfaceHopping

Independent electron surface hopping.

# References
- [Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009)](https://doi.org/10.1063/1.3125436)
- [Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009)](https://doi.org/10.1063/1.3122989)

"""
struct IESH_Tully{T} <: SurfaceHopping
    hopping_probability::Matrix{T}
    state::Vector{Int}
    new_state::Vector{Int}
    proposed_state::Vector{Int}
    n_electrons::Int
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    psi::Matrix{Complex{T}}

    function IESH_Tully{T}(states::Integer, n_electrons::Integer) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)
        psi = zeros(n_electrons*2+1,n_electrons)

        new{T}(hopping_probability, state, new_state, proposed_state,n_electrons, overlap, tmp, psi)
    end
end

function NQCDynamics.Simulation{IESH_Tully}(atoms::Atoms{S,T}, model::Model; n_electrons, kwargs...) where {S,T}
    NQCDynamics.Simulation(atoms, model, IESH_Tully{T}(NQCModels.nstates(model), n_electrons); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:IESH_Tully}, v, r)
    ψ = zeros(NQCModels.nstates(sim.calculator.model), sim.method.n_electrons)
    #occnum = zeros(NQCModels.nstates(sim.calculator.model), sim.method.n_electrons)
    phi = zeros(NQCModels.nstates(sim.calculator.model), sim.method.n_electrons)

    # Include Hamiltonian ground state into IESH propagation
    # should be initialized as eigenvector of the occupied states, which, for beginning are 
    # the lowest states n_electrons states
    # get eigenvectors
    Calculators.update_electronics!(sim.calculator, r)
    vecH = sim.calculator.eigen.vectors
    valH = sim.calculator.eigen.values
    state = collect(1:sim.method.n_electrons)
    ψ = vecH[:,state]
    phi = vecH[:,state]
    sim.method.psi .= vecH[:,state]
    
    # I now also define occnum, the occupation number. NOt sure if it will be necessary
    # for i=1:sim.method.n_electrons
    # #    Ψ[i,i-1] = 1
    #     occnum[i,i] = 1
    # end

    # This is in principle all right, but one needs to be careful with 
    # how the array is accessed.
    #state = collect(1:sim.method.n_electrons)

    SurfaceHoppingVariables(ComponentVector(v=v, r=r, σreal=ψ, σimag=zero(ψ)), state)
end

function DynamicsMethods.motion!(du, u, sim::Simulation{<:IESH_Tully}, t)
    dr = DynamicsUtils.get_positions(du)
    dv = DynamicsUtils.get_velocities(du)
    dσ = DynamicsUtils.get_quantum_subsystem(du)

    r = DynamicsUtils.get_positions(u)
    v = DynamicsUtils.get_velocities(u)
    σ = DynamicsUtils.get_quantum_subsystem(u)

    set_state!(u, sim.method.state) # Make sure the state variables match, 

    # Get velocities, I expect
    DynamicsUtils.velocity!(dr, v, r, sim, t)
    # Get nonadiabatic couplings and eigenvectors etc that correspond to present positions
    Calculators.update_electronics!(sim.calculator, r)
    # This only gets the new forces. It does not update the nuclear positions
    acceleration!(dv, v, r, sim, t, sim.method.state)
    #set_quantum_derivative!(sim.method.psi, v, σ, sim)
    # This is apparently needed, otherwise the program doesn't find it.
    propagate_wavefunction!(sim)
    #Internal way of propagation, needs correct wavefunction or density matrix
    println("psi")
    println(sim.method.psi)
    #set_quantum_derivative!(dσ, v, σ, sim)
    #println("sigma")
    #println(σ)
end

"""
Set the acceleration due to the force from the currently occupied states.
See Eq. 12 of Shenvi, Tully JCP 2009 paper.
"""
function acceleration!(dv, v, r, sim::Simulation{<:IESH_Tully}, t, state)
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
function set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:IESH_Tully})
    V = sim.calculator.eigen.values
    d = sim.calculator.nonadiabatic_coupling
    @views for i in axes(dσ, 2)       
        set_single_electron_derivative!(dσ[:,i], σ[:,i], V, v, d, sim.method.tmp)
    end
end
# Go to propagate the Schroedinger Equation.
# Needs to be the derivative, but I'm sabotaging the entire thing and just
# propagating the full favefunction
function propagate_wavefunction!(sim::Simulation{<:IESH_Tully})
    # Update the eigenvectors to the new positions before you propagate the Hamiltonian
    V = sim.calculator.eigen.values
    ab = zeros(Complex,size(sim.method.psi))
    
    Hvec = sim.calculator.eigen.vectors
    # Differences to Fortran code within numerical accuracy.
    ps = transpose(Hvec)*sim.method.psi

    # Set ps to value in Tully-IESH code to allow 1:1 comparison.
    # It checks out until here, though
    #ps = [0.99999999999999978+0.0im        1.00855884164624542E-016+0.0im; 1.00855884164624542E-016+0.0im       1.0000000000000000+0.0im;  -5.44007002986078857E-017+0.0im 2.23803759817048583E-016+0.0im;  1.85941304325904460E-016+0.0im      4.45375054458271422E-017+0.0im; 6.12028795037231673E-017+0.0im  8.48420904660896971E-017+0.0im]
    
    # The PES is originally in SI units with J, fs and Angstrom. 
    # In case the propagation is to be done in that unit system,
    # the following lines are needed. They give the same result.
    #hbar = 1.05457148e-34
    #dt = 1e-16 # timestep: 0.1 fs.
    #ac = exp.(-im*V*dt/hbar)
    #V1 = ustrip.(map( x -> auconvert(x*u"J"), V))

    # Timestep as 0.1 fs Hard-coded for now, since I have no idea how
    # to get it here.
    dt1 = ustrip(auconvert(1e-16*u"s"))
    ac = exp.(-im*V*dt1)
    
    for i in 1:sim.method.n_electrons
        for j in 1:(sim.method.n_electrons*2+1)
            ab[j,i] = ac[j]*ps[j,i]
        end
    end
    sim.method.psi .= Hvec*ab
    # Checks out until here.
end


"""
Hopping probability according to equation 21 in Shenvi, Roy, Tully 2009.
"""
function evaluate_hopping_probability!(sim::Simulation{<:IESH_Tully}, u, dt)

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    v = DynamicsUtils.get_velocities(u)

    # get phi_k (ϕk), the occupied orbitals
    vecH = sim.calculator.eigen.vectors
    phi = vecH[:,sim.method.state]
    psi = sim.method.psi
    #Get psi, the propagated Wavefunction of the entire system
    phipsi = det(transpose(phi)*psi)
    # Calculate Akk -- checks out against Tully code
    Akk = abs2(phipsi^2)

    #w This is wrong.
    #w S = sim.method.overlap
    proposed_state = sim.method.proposed_state
    prob = sim.method.hopping_probability
    # Get the nonadiabatic coupling.
    d = sim.calculator.nonadiabatic_coupling

    #w compute_overlap!(S, ψ, sim.method.state)

    #w det_current = det(S)
    #w Akk = abs2(det_current)
    pk= derivative(sim.calculator.model,DynamicsUtils.get_positions(u))
    #println(size(pk))
    #println(pk[:][1,1:2])
    
    fill!(prob, zero(eltype(prob)))
    for i in axes(prob, 2) # each electron
        for j in axes(prob, 1) # each basis state
            if j ∉ sim.method.state
                copyto!(proposed_state, sim.method.state)
                proposed_state[i] = j
                phij = vecH[:,proposed_state]
                # Akj (i.e. Aij) checks out mathematically against the fortran code
                # The numbers are quite different, but that seems to be a rounding
                # Problem again.
                Akj = calculate_Akj(sim.method.psi, phipsi, phij)
                # This looks correct. It adds for every atom.
                #Akj/Akk looks also reasonable and putting both into Real is, too, bcs AKK is a
                # The values for prob are quite small compared to the Fortran code.
                # This may either be numerical precision or a result of how the forces are treated
                # The formulae are correct.
                # For some reason, however, the forces have not changed, yet.
                for I in eachindex(v)
                    #prob[j,i] += 2 * v[I] * real(Akj/Akk) * d[I][sim.method.state[i], j] * dt
                    prob[j,i] += v[I] * d[I][sim.method.state[i], j]
                end
            end
        end
    end

    clamp!(prob, 0, 1) # Restrict probabilities between 0 and 1

    return nothing
end

"Equation 17 in Shenvi, Roy, Tully 2009. Uses equations 19 and 20."
function calculate_Akj(psi, phipsi, phij)
    
    ctemp = det(transpose(phij)*psi)
    # # Check Matrix-matrix multiplications
    # ab = zeros(Complex,(sim.method.n_electrons,sim.method.n_electrons))
    # for i in 1:sim.method.n_electrons
    #     for k in 1:(sim.method.n_electrons*2+1)
    #         for j in 1:sim.method.n_electrons
    #             ab[i,j] = ab[i,j]+transpose(phij)[i,k]*psi[k,j]
    #             #println(transpose(phij)[i,k]," ", psi[k,j]," ",  transpose(phij)[i,k]*psi[k,j]," ", ab[i,j]," ",i," ",j)
    #         end
    #     end
    # end
    return phipsi*conj(ctemp)
end

#w "Equation 20 in Shenvi, Roy, Tully 2009."
#w function compute_overlap!(S::Matrix, ψ, state)
#w     for i in axes(S, 2) # each electron
#w         for j in axes(S, 1) # each electron
#w             S[j,i] = ψ[state[j],i]
#w         end
#w     end
#w     return nothing
#w end

function select_new_state(sim::AbstractSimulation{<:IESH_Tully}, u)::Vector{Int}
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

function Estimators.diabatic_population(sim::Simulation{<:IESH_Tully}, u)
    Calculators.evaluate_potential!(sim.calculator, DynamicsUtils.get_positions(u))
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigen.vectors

    ψ = DynamicsUtils.get_quantum_subsystem(u)
    diabatic_ψ = zero(ψ)
    @views for i in axes(ψ, 2) # each electron
        diabatic_ψ[:,i] .= U*ψ[:,i]
    end
    diabatic_population = sum(abs2.(diabatic_ψ), dims=2)

    return diabatic_population
end

function Estimators.adiabatic_population(sim::Simulation{<:IESH_Tully}, u)
    population = zeros(NQCModels.nstates(sim.calculator.model))
    population[u.state] .= 1
    return population
end

function rescale_velocity!(sim::Simulation{<:IESH_Tully}, u)::Bool
    # ShakibHuo_JPhysChemLett_8_3073_2017_rescale
    new_state, old_state = symdiff(sim.method.new_state, sim.method.state)
    velocity = DynamicsUtils.get_velocities(u)
    
    # Calculate difference in eigenvalues between old and new state, weighed by mass:
    # c = (E_{new} - E_{old})/mass
    c = calculate_potential_energy_change(sim.calculator, new_state, old_state)
    # a = d^2_{new,old}/mass (where d is the nonadiabatic coupling)
    # b = v*d_{new,old}
    a, b = evaluate_a_and_b(sim, velocity, new_state, old_state)
    discriminant = b.^2 .- 2a.*c

    # If smaller zero, hop is frustrated
    any(discriminant .< 0) && return false

    root = sqrt.(discriminant)
    plus = (b .+ root) ./ a
    minus = (b .- root) ./ a 
    velocity_rescale = sum(abs.(plus)) < sum(abs.(minus)) ? plus : minus
    perform_rescaling!(sim, velocity, velocity_rescale, new_state, old_state)

    return true
end

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