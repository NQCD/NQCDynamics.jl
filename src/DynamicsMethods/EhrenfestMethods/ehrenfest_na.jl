using NQCModels: eachelectron, eachstate, mobileatoms, dofs
using LinearAlgebra: LinearAlgebra

struct EhrenfestNA{T}
    tmp::Vector{Complex{T}}
    quantum_propagator::Matrix{Complex{T}}
    tmp_matrix_complex_square1::Matrix{Complex{T}}
    tmp_matrix_complex_square2::Matrix{Complex{T}}
    tmp_matrix_complex_rect1::Matrix{Complex{T}}
    tmp_matrix_complex_rect2::Matrix{Complex{T}}
    function EhrenfestNA{T}(states::Integer, n_electrons::Integer) where {T}

        tmp = zeros(Complex{T}, states)
        quantum_propagator = zeros(Complex{T}, states, states)
        tmp_matrix_complex_square1 = zeros(Complex{T}, states, states)
        tmp_matrix_complex_square2 = zeros(Complex{T}, states, states)
        tmp_matrix_complex_rect1 = zeros(Complex{T}, states, n_electrons)
        tmp_matrix_complex_rect2 = zeros(Complex{T}, states, n_electrons)

        new{T}(tmp, quantum_propagator,
            tmp_matrix_complex_square1, tmp_matrix_complex_square2,
            tmp_matrix_complex_rect1, tmp_matrix_complex_rect2,
        )
    end
end

function NQCDynamics.Simulation{EhrenfestNA}(atoms::Atoms{T}, model::Model; kwargs...) where {T}
    NQCDynamics.Simulation(atoms, model, EhrenfestNA{T}(NQCModels.nstates(model), NQCModels.nelectrons(model)); kwargs...)
end

function DynamicsMethods.DynamicsVariables(sim::Simulation{<:EhrenfestNA}, v, r)
    ψ = zeros(NQCModels.nstates(sim), NQCModels.nelectrons(sim))
    
    for i in eachelectron(sim)
        ψ[i,i] = 1
    end

    return ComponentVector(v=v, r=r, σreal=ψ, σimag=zero(ψ))
end

function DynamicsUtils.acceleration!(dv, v, r, sim::Simulation{<:EhrenfestNA}, t, ψ)
    dv .= zero(eltype(dv))
    NQCModels.state_independent_derivative!(sim.calculator.model, dv, r)
    LinearAlgebra.lmul!(-1, dv)

    adiabatic_derivative = Calculators.get_adiabatic_derivative(sim.calculator, r)
    @inbounds for i in mobileatoms(sim)
        for j in dofs(sim)
            for electron in eachelectron(sim)
                for m in eachstate(sim)
                    # dv[j,i] -= adiabatic_derivative[j,i][m,m] * abs2(ψ[m,electron])
                    for n in eachstate(sim)
                        dv[j,i] -= adiabatic_derivative[j,i][n,m] * real(ψ[n,electron] * conj(ψ[m,electron]))
                    end
                end
            end
        end
    end
    DynamicsUtils.divide_by_mass!(dv, sim.atoms.masses)
    return nothing
end

function DynamicsUtils.set_quantum_derivative!(dσ, v, σ, sim::Simulation{<:EhrenfestNA})
    V = sim.calculator.eigen.values
    d = sim.calculator.nonadiabatic_coupling
    @views for i in eachelectron(sim)
        DynamicsUtils.set_single_electron_derivative!(dσ[:,i], σ[:,i], V, v, d, sim.method.tmp)
    end
end

function DynamicsUtils.classical_potential_energy(sim::Simulation{<:EhrenfestNA}, u)
    eigen = Calculators.get_eigen(sim.calculator, DynamicsUtils.get_positions(u))
    potential = NQCModels.state_independent_potential(sim.calculator.model, DynamicsUtils.get_positions(u))
    ψ = DynamicsUtils.get_quantum_subsystem(u)

    for electron in eachelectron(sim)
        for i in eachstate(sim)
            potential += eigen.values[i] * abs2(ψ[i,electron])
        end
    end
    return potential
end

function DynamicsUtils.classical_hamiltonian(sim::Simulation{<:EhrenfestNA}, u)
    kinetic = DynamicsUtils.classical_kinetic_energy(sim, DynamicsUtils.get_velocities(u))
    potential = DynamicsUtils.classical_potential_energy(sim, u)
    return kinetic + potential
end
