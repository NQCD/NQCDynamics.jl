export Langevin
export ThermalLangevin

struct Langevin{T<:AbstractFloat} <: Method
    γ::T
    σ::Matrix{T}
end

function Langevin{T}(γ, temperature, masses, DoFs) where {T}
    σ = sqrt.(temperature ./ repeat(Array(masses'), Int(DoFs), 1))
    Langevin(T(γ), T.(σ))
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::Simulation{<:Langevin})
    create_problem(u0.x, tspan, sim)
end
function create_problem(u0::ArrayPartition, tspan::Tuple, sim::Simulation{<:Langevin})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end
select_algorithm(sim::AbstractSimulation{<:Langevin}) = BAOAB(sim.method.γ)

function friction!(du, r, sim::AbstractSimulation{<:Langevin}, t)
    du .= sim.method.σ
end

struct ThermalLangevin{T<:Real} <: Method
    γ::T
end

function create_problem(u0::ClassicalDynamicals, tspan::Tuple, sim::AbstractSimulation{<:ThermalLangevin})
    DynamicalSDEProblem(acceleration!, velocity!, friction!, get_velocities(u0), get_positions(u0), tspan, sim)
end
select_algorithm(::RingPolymerSimulation{<:ThermalLangevin}) = BCOCB()
