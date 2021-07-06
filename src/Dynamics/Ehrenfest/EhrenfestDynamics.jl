"""
    $(TYPEDEF)

Abstract type for Ehrenfest method.

"""
abstract type AbstractEhrenfest <: Method end

function motion!(du, u, sim::AbstractSimulation{<:AbstractEhrenfest}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_quantum_subsystem(du)
    r = get_positions(u)
    v = get_velocities(u)
    σ = get_quantum_subsystem(u)
    velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    set_quantum_derivative!(dσ, v, σ, sim)
end

function set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:AbstractEhrenfest})
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.nonadiabatic_coupling[I]
    end
    mul!(sim.calculator.tmp_mat_complex1, V, σ)
    mul!(sim.calculator.tmp_mat_complex2, σ, V)
    @. dσ = -im * (sim.calculator.tmp_mat_complex1 - sim.calculator.tmp_mat_complex2)
end

function create_problem(u0, tspan, sim::AbstractSimulation{<:AbstractEhrenfest})
    ODEProblem(motion!, u0, tspan, sim)
end

include("ehrenfest_variables.jl")
include("ehrenfest.jl")
include("ehrenfest_rpmd.jl")