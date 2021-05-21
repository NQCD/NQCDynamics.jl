using StatsBase: mean
using .Calculators: DiabaticCalculator, RingPolymerDiabaticCalculator

export ehrenfest_basic

mutable struct ehrenfest_basic{T} <: Ehrenfest
    density_propagator::Matrix{Complex{T}}
    state::Int
    function ehrenfest_basic{T}(states::Integer) where {T}
        density_propagator = zeros(states, states)
        new{T}(density_propagator)
    end
end


function acceleration!(dv, v, r, sim::Simulation{<:ehrenfest_basic}, t, σ)
    for i in axes(dv, 2)
        for j in axes(dv, 1)
            dv[j,i] = round(-sum(sim.calculator.adiabatic_derivative[j,i] .* σ) / sim.atoms.masses[i])
        end
    end
    return nothing
end

"""
    get_population(sim::Simulation{<:ehrenfest_basic}, u)

"""
function get_population(sim::Simulation{<:ehrenfest_basic}, u)
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    U = sim.calculator.eigenvectors

    σ = copy(get_density_matrix(u))

    return real.(diag(U * σ * U'))
end

function NonadiabaticMolecularDynamics.evaluate_hamiltonian(sim::Simulation{<:ehrenfest_basic}, u)
    k = evaluate_kinetic_energy(sim.atoms.masses, get_velocities(u))
    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    Calculators.eigen!(sim.calculator)
    p = sim.calculator.eigenvalues[u.state]
    return k + p
end
