export NRPMD
export RingPolymerMappingVariables
export get_population

"""
$(TYPEDEF)

Nonadiabatic ring polymer molecular dynamics
"""
struct NRPMD{T} <: Method
    temp_q::Vector{T}
    temp_p::Vector{T}
    function NRPMD{T}(n_states::Integer) where {T}
        new{T}(zeros(n_states), zeros(n_states))
    end
end

"""
$(TYPEDEF)

Type for containing the classical variables for ring polymer mapping methods.
"""
mutable struct RingPolymerMappingVariables{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{RingPolymerArray{T},RingPolymerArray{T}, Matrix{T}, Matrix{T}}}
end

function RingPolymerMappingVariables(v::RingPolymerArray{T}, r::RingPolymerArray{T},
                                      n_states::Integer, state::Integer) where {T}

    n_beads = size(r)[3]
    qmap = zeros(T, n_states, n_beads)
    pmap = zeros(T, n_states, n_beads)

    θ = rand(T, n_states, n_beads) .* 2π
    qmap .= 1 ./ sqrt.(tan.(-θ) .^ 2 .+ 1)
    qmap[state,:] .*= sqrt(3)
    pmap = qmap .* tan.(-θ)

    RingPolymerMappingVariables{T}(ArrayPartition(v, r, pmap, qmap))
end

function RingPolymerMappingVariables(v::Matrix, r::Matrix, n_beads::Integer,
        n_states::Integer, state::Integer)
    v = RingPolymerArray(cat([v for i=1:n_beads]..., dims=3))
    r = RingPolymerArray(cat([r for i=1:n_beads]..., dims=3))
    RingPolymerMappingVariables(v, r, n_states, state)
end

function RingPolymerMappingVariables(v::Array{T,3}, r::Array{T,3}, n_states, state) where {T}
    RingPolymerMappingVariables(RingPolymerArray(v), RingPolymerArray(r), n_states, state)
end

get_mapping_positions(z::RingPolymerMappingVariables) = z.x.x[4]
get_mapping_momenta(z::RingPolymerMappingVariables) = z.x.x[3]
get_mapping_positions(z::RingPolymerMappingVariables, i::Integer) = @view z.x.x[4][:,i]
get_mapping_momenta(z::RingPolymerMappingVariables, i::Integer) = @view z.x.x[3][:,i]

function motion!(du::RingPolymerMappingVariables, u::RingPolymerMappingVariables,
        sim::RingPolymerSimulation{<:NRPMD}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u::RingPolymerMappingVariables,
        sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        for j in range(sim.atoms)
            for k=1:sim.DoFs
                D = sim.calculator.derivative[k,j,i]
                mul!(sim.method.temp_q, D, qmap)
                mul!(sim.method.temp_p, D, pmap)
                dv[k,j,i] = dot(qmap, sim.method.temp_q)
                dv[k,j,i] += dot(pmap, sim.method.temp_p)
                dv[k,j,i] -= tr(D)
                dv[k,j,i] /= -2sim.atoms.masses[j]
            end
        end
    end
    apply_interbead_coupling!(dv, get_positions(u), sim)
end

function set_mapping_force!(du::RingPolymerMappingVariables,
        u::RingPolymerMappingVariables, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        mul!(get_mapping_positions(du, i), V, get_mapping_momenta(u, i))
        mul!(get_mapping_momenta(du, i), V, get_mapping_positions(u, i))
        get_mapping_momenta(du, i) .*= -1
    end
end

function get_population(::RingPolymerSimulation{<:NRPMD}, u::RingPolymerMappingVariables)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    sum(qmap.^2 + pmap.^2 .- 1; dims=2) / 2size(qmap, 2)
end
