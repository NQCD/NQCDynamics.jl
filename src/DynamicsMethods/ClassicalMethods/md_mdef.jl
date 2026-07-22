#= using ACEglefriction

struct MD_MDEF{M,mm} <: AbstractMDEF
    mass_scaling::M
    memory_model::mm
    r_prime::Vector{Float64}
end

MD_MDEF(masses::AbstractVector, DoFs::Integer,memory_model) = MD_MDEF(get_mass_scale_matrix(masses, DoFs), memory_model)

function NQCDynamics.Simulation{MD_MDEF}(atoms::Atoms, model::ACEglefriction.GLEFrictionCalculator, memory_model; kwargs...)
    NQCDynamics.Simulation(atoms, model, MD_MDEF(atoms.masses, ndofs(model), memory_model); kwargs...)
end

function friction!(g, r, sim::AbstractSimulation{<:MD_MDEF}, t)
    position_to_vector!(sim.method.r_prime, r)
    r_prime = sim.method.r_prime
    g .= sim.calculator.model(r_prime)
end

function step_O!(integrator_cache::MD_MDEF_Cache, integrator)
    @unpack t, dt, W, p, sqdt = integrator
    @unpack dutmp, flatdutmp, tmp1, tmp2, gtmp, noise, half, c1, c2 = integrator_cache

    # gtmp = friction tensor
    # p = sim
    # tmp1 = momentum
    F = exp(-gtmp * dt)
    cd = cholesky(Symmetric(I - F * transpose(F)), RowMaximum(), check=true)
    S  = cd.L[invperm(cd.p), 1:cd.rank]
    m_sqrt = sqrt.(p.atoms.masses)
    # Work in mass-weighted coordinates z = (p/√m, s)

    tmp1 .= dutmp ./ p.atoms.masses
    ps = F * tmp1 ./ m_sqrt
    
end

function O_step(gc::GLEFrictionCalculator, h::T, x::Vector{T}, p::Vector{T}, s::Vector{T}; β::T=1.0) where {T}
    Γ = friction!(yadaydadaya)
    F = exp(-h * Γ)
    # Pivoted Cholesky of I - F Fᵀ; rank-revealing to handle numerical near-rank-deficiency
    cd = cholesky(Symmetric(I - F * transpose(F)), RowMaximum(), check=true)
    S  = cd.L[invperm(cd.p), 1:cd.rank]

    m_sqrt = sqrt.(get_masses(gc))
    # Work in mass-weighted coordinates z = (p/√m, s)
    ps = F * vcat(p ./ m_sqrt, s)
    # Add thermal noise in mass-weighted space
    ps += sqrt(1.0/β) .* (S * randn(T, size(S, 2)))
    # Restore physical momentum units for the p block
    ps[1:length(p)] = m_sqrt .* ps[1:length(p)]
    return ps[1:length(p)], ps[length(p)+1:end]
end =#