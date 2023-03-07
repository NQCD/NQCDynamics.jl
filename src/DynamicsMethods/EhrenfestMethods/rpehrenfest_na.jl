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
    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    @inbounds for b in axes(dv, 3) # eachbead
        for i in mobileatoms(sim) # eachatom
            for j in dofs(sim) # eachdof
                for electron in eachelectron(sim)
                    for m in eachstate(sim)
                        for n in eachstate(sim)
                            dv[j,i,b] -= adiabatic_derivative[j,i,b][n,m] * real(ψ[n,electron] * conj(ψ[m,electron]))
                        end
                    end
                end
            end
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)

    return nothing
end

function DynamicsUtils.classical_potential_energy(sim::RingPolymerSimulation{<:EhrenfestNA}, u)
    r = DynamicsUtils.get_positions(u)
    eigen = Calculators.get_eigen(sim.calculator, r)
    potential = NQCModels.state_independent_potential(sim.calculator.model, r)
    ψ = DynamicsUtils.get_quantum_subsystem(u)

    for b in axes(r, 3) # eachbead
        for electron in eachelectron(sim)
            for i in eachstate(sim)
                potential += eigen[b].values[i] * abs2(ψ[i,electron])
            end
        end
    end

    spring = RingPolymers.get_spring_energy(sim.beads, sim.atoms.masses, r)
    return potential + spring
end
