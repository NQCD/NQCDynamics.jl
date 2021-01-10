export NRPMD
export RingPolymerMappingPhasespace
export get_population

"""
    NRPMD <: Method

Nonadiabatic ring polymer molecular dynamics
"""
struct NRPMD <: Method end

"""
    RingPolymerMappingPhasespace{T} <: DynamicalVariables{T}

Type for containing the classical variables for ring polymer mapping methods.

Can be used for both NRPMD and MVRPMD.
"""
mutable struct RingPolymerMappingPhasespace{T} <: DynamicalVariables{T}
    x::ArrayPartition{T, Tuple{Array{T,3}, Array{T,3}, Matrix{T}, Matrix{T}}}
end

function RingPolymerMappingPhasespace(R::Array{T,3}, P::Array{T,3},
                                      n_states::Integer, state::Integer) where {T}

    n_beads = size(R)[3]
    qmap = zeros(T, n_states, n_beads)
    pmap = zeros(T, n_states, n_beads)

    θ = rand(T, n_states, n_beads) .* 2π
    qmap .= 1 ./ sqrt.(tan.(-θ) .^ 2 .+ 1)
    qmap[state,:] .*= sqrt(3)
    pmap = qmap .* tan.(-θ)

    RingPolymerMappingPhasespace{T}(ArrayPartition(R, P, qmap, pmap))
end

function RingPolymerMappingPhasespace(R::Matrix, P::Matrix, n_beads::Integer,
        n_states::Integer, state::Integer)
    R = cat([R for i=1:n_beads]..., dims=3)
    P = cat([P for i=1:n_beads]..., dims=3)
    RingPolymerMappingPhasespace(R, P, n_states, state)
end

get_mapping_positions(z::RingPolymerMappingPhasespace) = z.x.x[3]
get_mapping_momenta(z::RingPolymerMappingPhasespace) = z.x.x[4]
get_mapping_positions(z::RingPolymerMappingPhasespace, i::Integer) = @view z.x.x[3][:,i]
get_mapping_momenta(z::RingPolymerMappingPhasespace, i::Integer) = @view z.x.x[4][:,i]

function motion!(du::RingPolymerMappingPhasespace, u::RingPolymerMappingPhasespace,
        sim::RingPolymerSimulation{<:NRPMD}, t)
    set_velocity!(du, u, sim)
    set_force!(du, u, sim)
    apply_interbead_coupling!(du, u, sim)
    set_mapping_force!(du, u, sim)
end

function set_force!(du::RingPolymerMappingPhasespace, u::RingPolymerMappingPhasespace,
        sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_derivative!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        qmap = get_mapping_positions(u, i)
        pmap = get_mapping_momenta(u, i)
        for j in range(sim.atoms)
            for k=1:sim.DoFs
                D = sim.calculator.derivative[k,j,i]
                get_momenta(du)[k,j,i] = -(qmap'D*qmap + pmap'D*pmap - tr(D))/2
            end
        end
    end
end

function set_mapping_force!(du::RingPolymerMappingPhasespace,
        u::RingPolymerMappingPhasespace, sim::RingPolymerSimulation{<:NRPMD})

    Calculators.evaluate_potential!(sim.calculator, get_positions(u))
    for i in range(sim.beads)
        V = sim.calculator.potential[i]
        get_mapping_positions(du, i) .= V * get_mapping_momenta(u, i)
        get_mapping_momenta(du, i) .= -V * get_mapping_positions(u, i)
    end
end

function get_population(u::RingPolymerMappingPhasespace)
    qmap = get_mapping_positions(u)
    pmap = get_mapping_momenta(u)
    sum(qmap.^2 + pmap.^2 .- 1; dims=2) / 2size(qmap, 2)
end
