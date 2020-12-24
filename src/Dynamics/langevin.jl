export Langevin

struct Langevin{T<:AbstractFloat} <: Method
    η::T
end
Langevin(η::Integer) = Langevin{Float64}(η)

function set_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:Langevin})
    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    get_momenta(du) .= -sim.calculator.derivative .- sim.method.η .* get_momenta(u)
end

function random_force!(du::Phasespace, u::Phasespace, sim::Simulation{<:Langevin}, t)
    get_momenta(du) .= sqrt.(2sim.method.η * sim.temperature  .* sim.atoms.masses')
end
