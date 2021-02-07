export Langevin

struct Langevin{T<:AbstractFloat} <: Method
    η::T
end
Langevin(η::Integer) = Langevin{Float64}(η)

# function create_problem(u0::Phasespace, tspan::Tuple, sim::Simulation{<:Langevin})
#     SDEProblem(motion!, random_force!, u0, tspan, sim)
# end
# select_algorithm(::AbstractSimulation{<:Langevin}) = LambaEulerHeun()

# function set_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:Langevin})
#     Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
#     get_momenta(du) .= -sim.calculator.derivative .- sim.method.η .* get_momenta(u)
# end

# function random_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:Langevin}, t)
#     for i=1:sim.DoFs
#         get_momenta(du)[i,:] .= sqrt.(2sim.method.η * sim.temperature  .* sim.atoms.masses)
#     end
# end
