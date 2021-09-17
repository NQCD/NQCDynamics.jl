using NonadiabaticModels: nstates

using .DynamicsMethods.ClassicalMethods: Classical, Langevin, ThermalLangevin, MDEF
using .DynamicsMethods.MappingVariableMethods: NRPMD
using .DynamicsMethods.SurfaceHoppingMethods: FSSH, IESH
using .DynamicsMethods.EhrenfestMethods: Ehrenfest

function Simulation(atoms::Atoms, model::Model; kwargs...)
    Simulation{Classical}(atoms, model; kwargs...)
end
function RingPolymerSimulation(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    RingPolymerSimulation{Classical}(atoms, model, n_beads; kwargs...)
end

function Simulation{Classical}(atoms::Atoms, model::Model; kwargs...)
    Simulation(atoms, model, Classical(); kwargs...)
end
function RingPolymerSimulation{Classical}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...)
    RingPolymerSimulation(atoms, model, Classical(), n_beads::Integer; kwargs...)
end

function Simulation{MDEF}(atoms::Atoms, model::Model; DoFs=3, kwargs...)
    Simulation(atoms, model, MDEF(atoms.masses, DoFs); DoFs=DoFs, kwargs...)
end
function RingPolymerSimulation{MDEF}(atoms::Atoms, model::Model, n_beads::Integer; DoFs=3, kwargs...)
    RingPolymerSimulation(atoms, model, MDEF(atoms.masses, DoFs), n_beads; DoFs=DoFs, kwargs...)
end
    
function RingPolymerSimulation{NRPMD}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, NRPMD{T}(nstates(model)), n_beads; kwargs...)
end

function Simulation{Langevin}(atoms::Atoms{S,T}, model::Model; γ=1, temperature=0u"K", DoFs=3, kwargs...) where {S,T}
    Simulation(atoms, model, Langevin{T}(γ, austrip(temperature), atoms.masses, DoFs);
               temperature=temperature, DoFs=DoFs, kwargs...)
end
function RingPolymerSimulation{ThermalLangevin}(atoms::Atoms, model::Model, n_beads::Integer; γ=1, kwargs...)
    RingPolymerSimulation(atoms, model, ThermalLangevin(γ), n_beads; kwargs...)
end

function Simulation{FSSH}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T}
    Simulation(atoms, model, FSSH{T}(nstates(model)); kwargs...)
end
function RingPolymerSimulation{FSSH}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, FSSH{T}(nstates(model)), n_beads; kwargs...)
end

function Simulation{IESH}(atoms::Atoms{S,T}, model::Model; n_electrons, kwargs...) where {S,T}
    Simulation(atoms, model, IESH{T}(nstates(model), n_electrons); kwargs...)
end

function Simulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T}
    Simulation(atoms, model, Ehrenfest{T}(nstates(model)); kwargs...)
end
function RingPolymerSimulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}
    RingPolymerSimulation(atoms, model, Ehrenfest{T}(nstates(model)), n_beads; kwargs...)
end
