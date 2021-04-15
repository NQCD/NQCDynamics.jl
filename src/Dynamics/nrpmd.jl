export NRPMD
export RingPolymerMappingDynamicals
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

Can be used for both NRPMD and MVRPMD.
"""
mutable struct RingPolymerMappingDynamicals{T,D} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{D,D, Matrix{T}, Matrix{T}}}
end

function RingPolymerMappingDynamicals(v::RingPolymerArray{T}, r::RingPolymerArray{T},
                                      n_states::Integer, state::Integer) where {T}

    n_beads = size(r)[3]
    qmap = zeros(T, n_states, n_beads)
    pmap = zeros(T, n_states, n_beads)

    θ = rand(T, n_states, n_beads) .* 2π
    qmap .= 1 ./ sqrt.(tan.(-θ) .^ 2 .+ 1)
    qmap[state,:] .*= sqrt(3)
    pmap = qmap .* tan.(-θ)

    RingPolymerMappingDynamicals{T,RingPolymerArray{T}}(ArrayPartition(v, r, pmap, qmap))
end

function RingPolymerMappingDynamicals(v::Matrix, r::Matrix, n_beads::Integer,
        n_states::Integer, state::Integer)
    v = RingPolymerArray(cat([v for i=1:n_beads]..., dims=3))
    r = RingPolymerArray(cat([r for i=1:n_beads]..., dims=3))
    RingPolymerMappingDynamicals(v, r, n_states, state)
end

get_mapping_positions(z::RingPolymerMappingDynamicals) = z.x.x[4]
get_mapping_momenta(z::RingPolymerMappingDynamicals) = z.x.x[3]
get_mapping_positions(z::RingPolymerMappingDynamicals, i::Integer) = @view z.x.x[4][:,i]
get_mapping_momenta(z::RingPolymerMappingDynamicals, i::Integer) = @view z.x.x[3][:,i]

function motion!(du::RingPolymerMappingDynamicals, u::RingPolymerMappingDynamicals,
        sim::RingPolymerSimulation{<:NRPMD}, t)
    dr = get_positions(du)
    dv = get_velocities(du)
    r = get_positions(u)
    v = get_velocities(u)
    velocity!(dr, v, r, sim, t)
    acceleration!(dv, u, sim)
    set_mapping_force!(du, u, sim)
end

function acceleration!(dv, u::RingPolymerMappingDynamicals,
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

function set_mapping_force!(du::RingPolymerMappingDynamicals,
        u::RingPolymerMappingDynamicals, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        mul!(get_mapping_positions(du, i), V, get_mapping_momenta(u, i))
        mul!(get_mapping_momenta(du, i), V, get_mapping_positions(u, i))
        get_mapping_momenta(du, i) .*= -1
    end
end

function get_population(u::RingPolymerMappingDynamicals)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    sum(qmap.^2 + pmap.^2 .- 1; dims=2) / 2size(qmap, 2)
end
