using RingPolymerArrays: eachbead

function RingPolymerSimulation{AdiabaticIESH}(atoms::Atoms{T}, model::Model, n_beads::Integer;
    rescaling=:standard, estimate_probability=true, kwargs...
) where {T}
    RingPolymerSimulation(atoms, model,
        AdiabaticIESH{T}(
            NQCModels.nstates(model),
            NQCModels.nelectrons(model),
            rescaling,
            estimate_probability
        ),
        n_beads;
        kwargs...
    )
end

function DynamicsUtils.acceleration!(dv, v, r, sim::RingPolymerSimulation{<:AbstractIESH}, t, state)
    fill!(dv, zero(eltype(dv)))

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    @inbounds for b in axes(dv,3) 
        @views NQCModels.state_independent_derivative!(sim.calculator.model, dv[:,:,b], r[:,:,b])
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
    eigs = Calculators.get_eigen(sim.calculator, r)
    potential = zero(eltype(r))
    for b in axes(r,3) # eachbead
        potential += NQCModels.state_independent_potential(sim.calculator.model, view(r,:,:,b))
        for i in u.state
            potential += eigs[b].values[i]
        end
    end
    return potential
end
