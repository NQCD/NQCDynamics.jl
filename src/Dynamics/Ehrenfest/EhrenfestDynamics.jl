"""
    $(TYPEDEF)

Abstract type for all Ehrenfest methods.

Surface hopping methods follow the structure set out in this file.
The nuclear and electronic variables are propagated by the `motion!` function.
The surface hopping procedure is handled by the `HoppingCallback` which
uses the functions `check_hop!` and `execute_hop!` as its `condition` and `affect!`.

To add a new surface hopping scheme, you must create a new struct
and define methods for `evaluate_hopping_probability!`, `select_new_state`,
and `rescale_velocity!`.

See `fssh.jl` for an example implementation.
"""
abstract type Ehrenfest <: Method end

function motion!(du, u, sim::AbstractSimulation{<:Ehrenfest}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    dσ = get_density_matrix(du)
    r = get_positions(u)
    v = get_velocities(u)
    σ = get_density_matrix(u)
    velocity!(dr, v, r, sim, t)
    Calculators.update_electronics!(sim.calculator, r)
    acceleration!(dv, v, r, sim, t, σ)
    set_density_matrix_derivative!(dσ, v, σ, sim)
end

function set_density_matrix_derivative!(dσ, v, σ, sim::Simulation{<:Ehrenfest})
    V = sim.method.density_propagator

    V .= diagm(sim.calculator.eigenvalues)
    for I in eachindex(v)
        @. V -= im * v[I] * sim.calculator.nonadiabatic_coupling[I]
    end
    mul!(sim.calculator.tmp_mat_complex1, V, σ)
    mul!(sim.calculator.tmp_mat_complex2, σ, V)
    @. dσ = -im * (sim.calculator.tmp_mat_complex1 - sim.calculator.tmp_mat_complex2)
end

function create_problem(u0, tspan, sim::AbstractSimulation{<:Ehrenfest})
    sim.method.state = u0.state
    ODEProblem(motion!, u0, tspan, sim)
end

include("ehrenfest_variables.jl")
include("ehrenfest.jl")