
using LinearAlgebra: mul!, det

export DiabaticIESH

struct DiabaticIESH{T} <: AbstractIESH
    hopping_probability::Matrix{T}
    state::Vector{Int}
    new_state::Vector{Int}
    proposed_state::Vector{Int}
    n_electrons::Int
    overlap::Matrix{Complex{T}}
    tmp::Vector{Complex{T}}
    rescaling::Symbol
    function DiabaticIESH{T}(states::Integer, n_electrons::Integer, rescaling::Symbol) where {T}

        hopping_probability = zeros(Int, states, n_electrons)
        state = zeros(Int, n_electrons)
        new_state = zero(state)
        proposed_state = zero(state)
        tmp = zeros(Complex{T}, states)
        overlap = zeros(Complex{T}, n_electrons, n_electrons)

        new{T}(hopping_probability, state, new_state, proposed_state,n_electrons, overlap, tmp, rescaling)
    end
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:DiabaticIESH}, v, r)

    ψ = zeros(NonadiabaticModels.nstates(sim.calculator.model), sim.method.n_electrons)
    for i=1:sim.method.n_electrons
        ψ[i,i] = 1
    end
    state = collect(1:sim.method.n_electrons)

    Calculators.update_electronics!(sim.calculator, r)
    U = sim.calculator.eigen.vectors
    @views for i=1:sim.method.n_electrons
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

"Equation 20 in Shenvi, Roy, Tully 2009."
function compute_overlap!(sim::Simulation{<:DiabaticIESH}, S::Matrix, ψ, state)
    ϕ = @view sim.calculator.eigen.vectors[:,state]
    mul!(S, ϕ', ψ)
    return nothing
end

function Estimators.adiabatic_population(sim::Simulation{<:DiabaticIESH}, u)
    population = zeros(NonadiabaticModels.nstates(sim.calculator.model))
    Calculators.update_electronics!(sim.calculator, DynamicsUtils.get_positions(u))
    U = sim.calculator.eigen.vectors
    for i=1:sim.method.n_electrons
        adiabatic = U' * DynamicsUtils.get_quantum_subsystem(u)[:,i]
        population .+= abs2.(adiabatic)
    end
    return population
end
