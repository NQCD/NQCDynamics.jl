using RingPolymerArrays: eachbead

function RingPolymerSimulation{AdiabaticIESH}(atoms::Atoms{T}, model::Model, n_beads::Integer;
    rescaling=:standard, estimate_probability=true, disable_hopping=false, decoherence=DecoherenceCorrectionNone(), kwargs...
) where {T}
    RingPolymerSimulation(atoms, model,
        AdiabaticIESH{T}(
            NQCModels.nstates(model),
            NQCModels.nelectrons(model),
            rescaling,
            estimate_probability,
            disable_hopping,
            decoherence
        ),
        n_beads;
        kwargs...
    )
end

function DynamicsUtils.acceleration!(dv, v, r, sim::RingPolymerSimulation{<:AbstractIESH}, t, state)
    fill!(dv, zero(eltype(dv)))

    adiabatic_derivative = NQCCalculators.get_adiabatic_derivative(sim.cache, r)
    @inbounds for b in axes(dv,3) 
        @views NQCModels.state_independent_derivative!(sim.cache.model, dv[:,:,b], r[:,:,b])
        for i in mobileatoms(sim)
            for j in dofs(sim)
                for k in state
                    # Contribution to the force from each occupied state `k`
                    dv[j,i,b] += adiabatic_derivative[j,i,b][k, k]
                end
            end
        end
    end
    lmul!(-1, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:AbstractIESH}, u)
    r = DynamicsUtils.get_positions(u)
    eigs = NQCCalculators.get_eigen(sim.cache, r)
    potential = zero(eltype(r))
    for b in axes(r,3) # eachbead
        potential += NQCModels.state_independent_potential(sim.cache.model, view(r,:,:,b))
        for i in u.state
            potential += eigs[b].values[i]
        end
    end
    return potential
end
