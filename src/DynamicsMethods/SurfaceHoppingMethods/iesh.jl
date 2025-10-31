using LinearAlgebra: mul!, det, Hermitian, diagm
using Unitful, UnitfulAtomic
using MuladdMacro: @muladd
using NQCDynamics: FastDeterminant
using NQCModels: eachelectron, eachstate, mobileatoms, dofs
using NQCDistributions: FermiDiracState, Adiabatic, Diabatic
using StatsBase: sample, Weights
using FastLapackInterface: FastLapackInterface
using SciMLBase: SciMLBase
using RecursiveArrayTools
using SciMLOperators

export AdiabaticIESH

abstract type AbstractIESH <: SurfaceHopping end

"""
    IESH{T} <: SurfaceHopping

Independent electron surface hopping.

# References
- [Shenvi, Roy, Tully, J. Chem. Phys. 130, 174107 (2009)](https://doi.org/10.1063/1.3125436)
- [Roy, Shenvi, Tully, J. Chem. Phys. 130, 174716 (2009)](https://doi.org/10.1063/1.3122989)

"""
struct AdiabaticIESH{T,D} <: AbstractIESH
    hopping_probability::Matrix{T}
    state::Vector{Int}
    new_state::Vector{Int}
    proposed_state::Vector{Int}
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    rescaling::Symbol
    quantum_propagator::Matrix{Complex{T}}
    LUws::FastLapackInterface.LUWs
    v_dot_d::Matrix{T}
    unoccupied::Vector{Int}
    estimate_probability::Bool
    disable_hopping::Bool
    decoherence::D
    function AdiabaticIESH{T}(states::Integer, n_electrons::Integer, rescaling::Symbol, estimate_probability::Bool, disable_hopping::Bool, decoherence) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)
        quantum_propagator = zeros(Complex{T}, states, states)
        LUws = FastLapackInterface.LUWs(n_electrons)
        v_dot_d = zeros(T, states, n_electrons)
        unoccupied = zeros(Int, states - n_electrons)

        new{T,typeof(decoherence)}(hopping_probability, state, new_state, proposed_state, overlap, tmp, rescaling, quantum_propagator,
            LUws, v_dot_d, unoccupied, estimate_probability, disable_hopping
        )
    end
end

unoccupied_states(sim::AbstractSimulation{<:AbstractIESH}) = sim.method.unoccupied

function NQCDynamics.Simulation{IESH_type}(atoms::Atoms{T}, model::Model;
    rescaling=:standard, estimate_probability=true, disable_hopping=false, decoherence=DecoherenceCorrectionNone(),
    kwargs...) where {T,IESH_type<:AbstractIESH}
    NQCDynamics.Simulation(
        atoms, model,
        IESH_type{T}(
            NQCModels.nstates(model),
            NQCModels.nelectrons(model),
            rescaling,
            estimate_probability,
            disable_hopping,
            decoherence
        ); 
        kwargs...
    )
end

function DynamicsMethods.DynamicsVariables(sim::AbstractSimulation{<:AdiabaticIESH}, v, r)
    ψ = zeros(NQCModels.nstates(sim), NQCModels.nelectrons(sim))
    
    for i in eachelectron(sim)
        ψ[i,i] = 1
    end
    state = collect(eachelectron(sim))
    SurfaceHoppingVariables(v=Float64.(v), r=Float64.(r), σreal=Float64.(ψ), σimag=Float64.(zero(ψ)), state=Float64.(state))
end

function DynamicsMethods.DynamicsVariables(sim::AbstractSimulation{<:AdiabaticIESH}, v, r, electronic::FermiDiracState{Adiabatic})
    tmp_v = copy(v)
    tmp_r = copy(r)
    
    ef_model = NQCModels.fermilevel(sim)
    ef_distribution = electronic.fermi_level
    ef_model ≈ ef_distribution || throw(error(
        """
        Fermi level of model and distribution do not match:
            Distribution: $(ef_distribution)
            Model: $(ef_model)
        Change one of them to make them the same.
        """
    ))

    # NQCDynamics.NQCCalculators.update_cache!(sim.cache, tmp_r) # Ensure hopping eigenvalues are up to date
    eigenvalues = DynamicsUtils.get_hopping_eigenvalues(sim, tmp_r)


    available_states = DynamicsUtils.get_available_states(electronic.available_states, NQCModels.nstates(sim))
    state = DynamicsUtils.sample_fermi_dirac_distribution(eigenvalues, NQCModels.nelectrons(sim), available_states, electronic.β)
    # state = DynamicsUtils.sample_fermi_dirac_distribution(eigenvalues, NQCModels.nelectrons(sim), available_states, electronic.β, electronic.fermi_level)

    ψ = zeros(NQCModels.nstates(sim), NQCModels.nelectrons(sim))
    for (i, j) in enumerate(state)
        ψ[j,i] = 1
    end

    SurfaceHoppingVariables(v=Float64.(tmp_v), r=Float64.(tmp_r), σreal=Float64.(ψ), σimag=Float64.(zero(ψ)), state=Float64.(state))
end

function DynamicsMethods.create_problem(u0, tspan, sim::AbstractSimulation{<:AbstractIESH})
    set_state!(sim.method, convert(Vector{Int},u0.state), sim) # state 
    set_unoccupied_states!(sim)

    electronic = SurfaceHoppingVariables(σreal = u0.σreal, σimag = u0.σimag, state = u0.state)
    nuclear = NamedArrayPartition(v = u0.v, r=u0.r)

    prob1 = OrdinaryDiffEq.DynamicalODEProblem(DynamicsUtils.acceleration!, DynamicsUtils.velocity!, u0.v, u0.r, tspan, (sim, deepcopy(electronic)))
    A = SciMLOperators.MatrixOperator(complex(ones(size(sim.cache.potential))), update_func! = DynamicsUtils.get_quantum_propagator!)
    prob2 = OrdinaryDiffEq.ODEProblem(A, electronic, tspan, (sim, deepcopy(nuclear)))

    return DynamicsMethods.IntegrationAlgorithms.CoupledODEProblem(prob1, prob2, DynamicsMethods.get_callbacks(sim))
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:AdiabaticIESH}, v, r, electronic::FermiDiracState{Diabatic})
    ef_model = NQCModels.fermilevel(sim)
    ef_distribution = electronic.fermi_level
    ef_model ≈ ef_distribution || throw(error(
        """
        Fermi level of model and distribution do not match:
            Distribution: $(ef_distribution)
            Model: $(ef_model)
        Change one of them to make them the same.
        """
    ))

    #NQCDynamics.NQCCalculators.update_cache!(sim.cache, r) # Ensure potential is up to date
    available_states = get_available_states(electronic.available_states, NQCModels.nstates(sim))

    potential = NQCCalculators.get_potential(sim.cache, r)
    energies = @view potential[diagind(potential)]

    available_states = get_available_states(electronic.available_states, NQCModels.nstates(sim))
    diabatic_state = sample_fermi_dirac_distribution(energies, NQCModels.nelectrons(sim), available_states, electronic.β)

    U = DynamicsUtils.evaluate_transformation(sim.cache, r)

    ψ = zeros(NQCModels.nstates(sim), NQCModels.nelectrons(sim))
    adiabatic_state = zeros(Int, NQCModels.nelectrons(sim))

    for i in eachindex(adiabatic_state)
        diabatic_population = zeros(NQCModels.nstates(sim))
        diabatic_population[diabatic_state[i]] = 1
        adiabatic_density = U' * diagm(diabatic_population) * U
        adiabatic_population = diag(adiabatic_density)

        while true
            s = sample(Weights(adiabatic_population))
            if !(s in adiabatic_state)
                adiabatic_state[i] = s
                break
            end
        end

        adiabatic_coefficients = adiabatic_density[:,1] ./ sqrt(adiabatic_population[1])
        ψ[:,i] .= adiabatic_coefficients
    end

    sort!(adiabatic_state)
    SurfaceHoppingVariables(v=Float64.(v), r=Float64.(r), σreal=Float64.(ψ), σimag=Float64.(zero(ψ)), state=Float64.(adiabatic_state))
end

"""
Set the acceleration due to the force from the currently occupied states.
See Eq. 12 of Shenvi, Tully JCP 2009 paper.
"""
function DynamicsUtils.acceleration!(dv, v, r, parameters::Tuple{Simulation{<:AbstractIESH}, Any}, t)
    sim, electronic = parameters

    NQCCalculators.NQCCalculators.update_cache!(sim.cache, r)
    state = sim.method.state
    dv .= zero(eltype(dv))
    NQCModels.state_independent_derivative!(sim.cache.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)

    #NQCDynamics.NQCCalculators.update_cache!(sim.cache, r)
    adiabatic_derivative = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
    @inbounds for i in mobileatoms(sim)
        for j in dofs(sim)
            for k in state
                # Contribution to the force from each occupied state `k`
                dv[j,i] -= adiabatic_derivative[j,i][k, k]
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
function DynamicsUtils.set_quantum_derivative!(du, u, parameters::Tuple{AbstractSimulation{<:AdiabaticIESH},Any}, t)

    sim, nuclei = parameters
    r = DynamicsUtils.get_positions(nuclei)
    v = DynamicsUtils.get_velocities(nuclei)
    
    NQCCalculators.update_cache!(sim.cache, r)
    σ = DynamicsUtils.get_quantum_subsystem(u)
    dσ = DynamicsUtils.get_quantum_subsystem(du)
    V = DynamicsUtils.get_hopping_eigenvalues(sim, r)
    d = DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r)
    @views for i in eachelectron(sim)
        DynamicsUtils.set_single_electron_derivative!(dσ[:,i], σ[:,i], V, v, d, sim.method.tmp)
    end
end
"""
Hopping probability according to equation 21 in Shenvi, Roy, Tully 2009.
"""
function evaluate_hopping_probability!(sim::AbstractSimulation{<:AbstractIESH}, u, nuclei, dt, random)
    ψ = DynamicsUtils.get_quantum_subsystem(u)
    v = DynamicsUtils.get_velocities(nuclei)

    S = sim.method.overlap
    proposed_state = sim.method.proposed_state
    prob = sim.method.hopping_probability
    r = DynamicsUtils.get_positions(nuclei)
    d = DynamicsUtils.get_hopping_nonadiabatic_coupling(sim, r)

    compute_overlap!(sim, S, ψ, sim.method.state)

    det_current = FastDeterminant.det!(S, sim.method.LUws)
    Akk = abs2(det_current)
    prefactor = 2dt / Akk

    fill!(prob, zero(eltype(prob)))

    evaluate_v_dot_d!(sim, v, d)

    if sim.method.estimate_probability
        estimate = prefactor * sum(abs, sim.method.v_dot_d) * (abs(real(det_current)) + abs(imag(det_current)))
        estimate < random && return
    end

    @inbounds for n in eachelectron(sim)
        for m in unoccupied_states(sim)
            copy!(proposed_state, sim.method.state)
            proposed_state[n] = m

            Akj = calculate_Akj(sim, S, ψ, det_current, proposed_state)

            probability = prefactor * real(Akj) * sim.method.v_dot_d[m,n]
            prob[m,n] = probability
        end
    end

    clamp!(prob, 0, 1) # Restrict probabilities between 0 and 1

    @debug begin
        maximum_probability = sum(prob)
        maximum_probability > 1 && @warn "Hopping probability is large, consider reducing the time step" maximum_probability

        if sim.method.estimate_probability
            if !(maximum_probability < estimate)
                @warn "Hopping probability exceeded estimate!" maximum_probability estimate det_current Akk
                error("The hopping probability should never exceed the estimate.")
            end
        end
    end

    return nothing
end

function evaluate_v_dot_d!(sim::AbstractSimulation{<:AbstractIESH}, v, d)
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

"Equation 17 in Shenvi, Roy, Tully 2009. Uses equations 19 and 20."
function calculate_Akj(sim::AbstractSimulation{<:AbstractIESH}, S::AbstractMatrix, ψ::AbstractMatrix, detS::Number, new_state::Vector)
    compute_overlap!(sim, S, ψ, new_state)
    det_new = FastDeterminant.det!(S, sim.method.LUws)

    return detS*conj(det_new)
end

"Equation 20 in Shenvi, Roy, Tully 2009."
function compute_overlap!(sim::AbstractSimulation{<:AdiabaticIESH}, S::Matrix, ψ, state)
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
    ψ = DynamicsUtils.get_quantum_subsystem(u).re

    eigen = NQCCalculators.get_eigen(sim.cache, DynamicsUtils.get_positions(u))
    transformation = eigen.vectors

    return iesh_diabatic_population(ψ, transformation, sim.method.state)
end

"""
Calculate the diabatic population in J. Chem. Theory Comput. 2022, 18, 4615−4626.
Eqs. 12, 13 and 14 describe the steps to calculat it though the code is written to use matrix operations
instead of summations to calculate the result.
"""
function iesh_diabatic_population(ψ::AbstractMatrix, transformation::AbstractMatrix, discrete_occupations::AbstractVector)
    nstates = size(ψ, 1)
    nelectrons = size(ψ, 2)

    population = zeros(nstates)
    single_electron_density_matrix = zeros(nstates, nstates)
    for i in 1:nelectrons
        c = view(ψ, :, i) # Single electron wavefunction coefficients

        mul!(single_electron_density_matrix, c, c')
        single_electron_density_matrix[diagind(single_electron_density_matrix)] .= 0

        occupied_state = discrete_occupations[i]
        single_electron_density_matrix[occupied_state, occupied_state] = 1

        population += diag(transformation * single_electron_density_matrix * transformation')
    end
    return population
end

function Estimators.adiabatic_population(sim::Simulation{<:AdiabaticIESH}, u)
    population = zeros(NQCModels.nstates(sim.cache.model))
    population[round.(Int, u.state)] .= 1
    return population
end

unpack_states(sim::AbstractSimulation{<:AbstractIESH}) = symdiff(sim.method.new_state, sim.method.state)#

ishoppingdisabled(method::AbstractIESH) = method.disable_hopping

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:AbstractIESH}, u)
    NQCCalculators.update_cache!(sim.cache, DynamicsUtils.get_positions(u)) # Ensure eigen is populated
    eigen = NQCCalculators.get_eigen(sim.cache, DynamicsUtils.get_positions(u))
    potential = NQCModels.state_independent_potential(sim.cache.model, DynamicsUtils.get_positions(u))
    for i in sim.method.state
        potential += eigen.values[i]
    end
    return potential
end

function iesh_check_hop!(nuclei, electrons, t, integrator)::Bool
    sim = integrator.int2.p[1]

    ishoppingdisabled(sim.method) && return false
    random = rand()
    evaluate_hopping_probability!(sim, electrons, nuclei, OrdinaryDiffEq.get_proposed_dt(integrator.int2), random)
    set_new_state!(sim.method, select_new_state(sim, electrons, random))
    bool = sim.method.new_state != sim.method.state
    return bool
end

iesh_check_hop!(u, t, integrator) = false

function iesh_execute_hop!(nuclei, electrons, integrator)
    sim = integrator.int2.p[1]
    if rescale_velocity!(sim, electrons, nuclei)
        set_state!(electrons, sim.method.new_state, sim)
        set_state!(sim.method, sim.method.new_state, sim)
        set_unoccupied_states!(sim)
    end
    return nothing
end

function set_unoccupied_states!(sim::AbstractSimulation{<:AbstractIESH})
    DynamicsUtils.set_unoccupied_states!(sim.method.unoccupied, sim.method.state)
end

const IESHCallback = DiffEqBase.DiscreteCallback(iesh_check_hop!, iesh_execute_hop!;
                                                    save_positions=(false, false))

function DynamicsMethods.get_callbacks(sim::AbstractSimulation{<:AbstractIESH})
    if sim.method.decoherence isa DecoherenceCorrectionEDC
        return SciMLBase.CallbackSet(IESHCallback, IESHDecoherenceEDC())
    else
        return SciMLBase.CallbackSet(IESHCallback)
    end
end

function IESHDecoherenceEDC()
    return DiffEqBase.DiscreteCallback(
        (u,t,integrator)->true,
        iesh_apply_decoherence_correction_edc!;
        save_positions=(false, false)
    )
end

function iesh_apply_decoherence_correction_edc!(integrator)
    sim = integrator.p
    dt = OrdinaryDiffEq.get_proposed_dt(integrator)
    eigen = NQCCalculators.get_eigen(sim.cache, DynamicsUtils.get_positions(integrator.u))
    Ekin = DynamicsUtils.classical_kinetic_energy(sim, integrator.u)
    ψ = DynamicsUtils.get_quantum_subsystem(integrator.u)
    @views for (i, state) in enumerate(sim.method.state)
        apply_decoherence_correction!(ψ[:,i], sim.method.decoherence, state, dt, eigen.values, Ekin)
    end
end

function set_state!(container, new_state::AbstractVector, sim::AbstractSimulation{<:AbstractIESH})
    container.state .= new_state
end