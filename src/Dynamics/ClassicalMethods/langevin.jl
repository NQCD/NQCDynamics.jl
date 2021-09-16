
struct Langevin{T<:AbstractFloat} <: Dynamics.Method
    γ::T
    σ::Matrix{T}
end

function Langevin{T}(γ, temperature, masses, DoFs) where {T}
    σ = sqrt.(temperature ./ repeat(Array(masses'), Int(DoFs), 1))
    Langevin(T(γ), T.(σ))
end

function create_problem(u0, tspan::Tuple, sim::Simulation{<:Langevin})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end
select_algorithm(sim::AbstractSimulation{<:Langevin}) = StochasticDiffEq.BAOAB(sim.method.γ)

function friction!(du, r, sim::AbstractSimulation{<:Langevin}, t)
    du .= sim.method.σ
end

struct ThermalLangevin{T<:Real} <: Dynamics.Method
    γ::T
end

function DynamicsVariables(::AbstractSimulation{<:Union{ThermalLangevin, Langevin}}, v, r)
    ArrayPartition(v, r)
end

function create_problem(u0, tspan::Tuple, sim::AbstractSimulation{<:ThermalLangevin})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end
Dynamics.select_algorithm(::RingPolymerSimulation{<:ThermalLangevin}) = Algorithms.BCOCB()
