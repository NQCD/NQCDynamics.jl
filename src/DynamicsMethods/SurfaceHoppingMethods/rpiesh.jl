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
    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    @inbounds for b in axes(dv,3) 
        for i in mobileatoms(sim)
            for j in dofs(sim)
                for k in state
                    # Contribution to the force from each occupied state `k`
                    dv[j,i,b] -= adiabatic_derivative[j,i,b][k, k]
                end
            end
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end
