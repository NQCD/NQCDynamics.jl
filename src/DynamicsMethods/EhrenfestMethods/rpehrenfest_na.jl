using NQCDynamics: RingPolymers

function RingPolymerSimulation{EhrenfestNA}(atoms::Atoms{T}, model::Model, n_beads::Integer; kwargs...) where {T}
    return RingPolymerSimulation(
        atoms,
        model,
        EhrenfestNA{T}(NQCModels.nstates(model), NQCModels.nelectrons(model)),
        n_beads;
        kwargs...
    )
end

function DynamicsUtils.acceleration!(dv, v, r, sim::RingPolymerSimulation{<:EhrenfestNA}, t, ψ)
    fill!(dv, zero(eltype(dv)))

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    @inbounds for b in axes(dv, 3) # eachbead
        @views NQCModels.state_independent_derivative!(sim.calculator.model, dv[:,:,b], r[:,:,b])
        for i in mobileatoms(sim) # eachatom
            for j in dofs(sim) # eachdof
                for electron in eachelectron(sim)
                    for m in eachstate(sim)
                        for n in eachstate(sim)
                            dv[j,i,b] += adiabatic_derivative[j,i,b][n,m] * real(ψ[n,electron] * conj(ψ[m,electron]))
                        end
                    end
                end
            end
        end
    end
    LinearAlgebra.lmul!(-1, dv)
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)

    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:EhrenfestNA}, u)
    r = DynamicsUtils.get_positions(u)
    eigen = Calculators.get_eigen(sim.calculator, r)
    ψ = DynamicsUtils.get_quantum_subsystem(u)

    potential = zero(eltype(r))
    for b in axes(r, 3) # eachbead
        potential += NQCModels.state_independent_potential(sim.calculator.model, view(r,:,:,b))
        for electron in eachelectron(sim)
            for i in eachstate(sim)
                potential += eigen[b].values[i] * abs2(ψ[i,electron])
            end
        end
    end
    return potential
end
