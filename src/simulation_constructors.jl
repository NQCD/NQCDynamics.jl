
Simulation(atoms::Atoms, model::Model; kwargs...) =
    Simulation{Classical}(atoms, model; kwargs...)

RingPolymerSimulation(atoms::Atoms, model::Model, n_beads::Integer; kwargs...) =
    RingPolymerSimulation{Classical}(atoms, model, n_beads; kwargs...)

Simulation{Classical}(atoms::Atoms, model::Model; kwargs...) =
    Simulation(atoms, model, Classical(); kwargs...)

RingPolymerSimulation{Classical}(atoms::Atoms, model::Model, n_beads::Integer; kwargs...) =
    RingPolymerSimulation(atoms, model, Classical(), n_beads::Integer; kwargs...)

Simulation{MDEF}(atoms::Atoms, model::Model; DoFs=3, kwargs...) =
    Simulation(atoms, model, MDEF(atoms.masses, DoFs); DoFs=DoFs, kwargs...)

RingPolymerSimulation{MDEF}(atoms::Atoms, model::Model, n_beads::Integer; DoFs=3, kwargs...) =
    RingPolymerSimulation(atoms, model, MDEF(atoms.masses, DoFs), n_beads; DoFs=DoFs, kwargs...)
    
RingPolymerSimulation{NRPMD}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T}=
    RingPolymerSimulation(atoms, model, NRPMD{T}(model.n_states), n_beads; kwargs...)

Simulation{CMM2}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T} = 
    Simulation(atoms, model, CMM2{T}(model.n_states); kwargs...)

Simulation{Langevin}(atoms::Atoms{S,T}, model::Model; γ=1, temperature=0u"K", DoFs=3, kwargs...) where {S,T} =
    Simulation(atoms, model, Langevin{T}(γ, austrip(temperature), atoms.masses, DoFs);
               temperature=temperature, DoFs=DoFs, kwargs...)

# FSSH is defined in src/Dynamics/SurfaceHopping
Simulation{FSSH}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T} =
    Simulation(atoms, model, FSSH{T}(model.n_states); kwargs...)

# IESH is defined in src/Dynamics/SurfaceHopping
Simulation{IESH}(atoms::Atoms{S,T}, model::Model; n_electrons, kwargs...) where {S,T} =
    Simulation(atoms, model, IESH{T}(model.n_states, n_electrons); kwargs...)

RingPolymerSimulation{FSSH}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T} =
    RingPolymerSimulation(atoms, model, FSSH{T}(model.n_states), n_beads; kwargs...)

RingPolymerSimulation{ThermalLangevin}(atoms::Atoms, model::Model, n_beads::Integer; γ=1, kwargs...) =
    RingPolymerSimulation(atoms, model, ThermalLangevin(γ), n_beads; kwargs...)

Simulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model; kwargs...) where {S,T} =
    Simulation(atoms, model, Ehrenfest{T}(model.n_states); kwargs...)

RingPolymerSimulation{Ehrenfest}(atoms::Atoms{S,T}, model::Model, n_beads::Integer; kwargs...) where {S,T} =
    RingPolymerSimulation(atoms, model, Ehrenfest{T}(model.n_states), n_beads; kwargs...)
