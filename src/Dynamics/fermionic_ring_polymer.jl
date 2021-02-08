export FermionicBath

struct FermionicBath{T<:AbstractFloat} <: Method
    η::T
    kF::T
    function FermionicBath(U::T, ρ::T, kF::T) where {T}
        η = 4/3 * π * U^2 * ρ * kF
        new{T}(η, kF)
    end
end

function NonadiabaticMolecularDynamics.evaluate_potential_energy(sim::RingPolymerSimulation{<:FermionicBath}, R::Array{T, 3}) where {T}
    Calculators.evaluate_potential!(sim.calculator, R)
    get_spring_energy(sim, R) + sum(sim.calculator.potential)[1] + get_bath_energy(sim, R)
end

function get_bath_energy(sim::RingPolymerSimulation{<:FermionicBath}, R::Array{T,3}) where {T}
    E = 0.0
    @views for i=1:length(sim.beads)
        for j=i+1:length(sim.beads)
            for k in sim.beads.quantum_atoms # Only for quantum nuclei
                E -= sinc(sim.method.kF*norm(R[:,k,i].-R[:,k,j])/π)^2 / dist(length(sim.beads), i, j)^2
            end
        end
    end
    E * 3/2 * sim.method.η * sim.temperature * length(sim.beads) / π / sim.method.kF^2
end

dist(n, i, j) = min(abs(i-j), n - abs(i-j))