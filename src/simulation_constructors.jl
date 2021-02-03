
Simulation(atoms::Atoms, model::Model; kwargs...) =
    Simulation{Classical}(atoms, model; kwargs...)

Simulation{Classical}(atoms::Atoms, model::Model; kwargs...) =
    Simulation(atoms, model, Classical(); kwargs...)

RingPolymerSimulation{Classical}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...) =
    RingPolymerSimulation(atoms, model, Classical(), n_beads::Integer; kwargs...)

Simulation{MDEF}(atoms::Atoms{S,T}, model::Model; DoFs=3, kwargs...) where {S,T} =
    Simulation(atoms, model, MDEF{T}(length(atoms), DoF=DoFs); DoFs=DoFs, kwargs...)

RingPolymerSimulation{NRPMD}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}=
    RingPolymerSimulation(atoms, model, NRPMD{T}(model.n_states), n_beads; kwargs...)

Simulation{Langevin}(atoms::Atoms{S,T}, model::Model; η=1, kwargs...) where {S,T} =
    Simulation(atoms, model, Langevin{T}(η); kwargs...)

Simulation{FSSH}(atoms::Atoms{S,T}, model::Model; DoFs=3, kwargs...) where {S,T} =
    Simulation(atoms, model, FSSH{T}(DoFs, length(atoms), model.n_states); DoFs=DoFs, kwargs...)