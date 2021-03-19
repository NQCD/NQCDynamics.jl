
Simulation(atoms::Atoms, model::Model; kwargs...) =
    Simulation{Classical}(atoms, model; kwargs...)

Simulation{Classical}(atoms::Atoms, model::Model; kwargs...) =
    Simulation(atoms, model, Classical(); kwargs...)

RingPolymerSimulation{Classical}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...) =
    RingPolymerSimulation(atoms, model, Classical(), n_beads::Integer; kwargs...)

Simulation{MDEF}(atoms::Atoms, model::Model; kwargs...) =
    Simulation(atoms, model, MDEF(); kwargs...)
    
Simulation{TwoTemperatureMDEF}(atoms::Atoms, model::Model, T::Function; kwargs...) =
    Simulation(atoms, model, TwoTemperatureMDEF(T); kwargs...)

RingPolymerSimulation{NRPMD}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}=
    RingPolymerSimulation(atoms, model, NRPMD{T}(model.n_states), n_beads; kwargs...)

Simulation{Langevin}(atoms::Atoms{S,T}, model::Model; γ=1, temperature=0u"K", DoFs=3, kwargs...) where {S,T} =
    Simulation(atoms, model, Langevin{T}(γ, austrip(temperature), atoms.masses, DoFs);
               temperature=temperature, DoFs=DoFs, kwargs...)

Simulation{FSSH}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T} =
    Simulation(atoms, model, FSSH{T}(length(atoms), model.n_states); kwargs...)

RingPolymerSimulation{FSSH}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T} =
    RingPolymerSimulation(atoms, model, FSSH{T}(length(atoms), model.n_states), n_beads; kwargs...)
