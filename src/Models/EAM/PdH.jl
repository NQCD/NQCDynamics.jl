export PdH

using LinearAlgebra
using Unitful
using FiniteDiff

"""
Ciufo, R.A., Henkelman, G. Embedded atom method potential for hydrogen on palladium surfaces. J Mol Model 26, 336 (2020).
"""
struct PdH <: AdiabaticModel
    χ::Float64
    Φ::Float64
    δ::Float64
    β_PdPd::Float64
    η::Float64
    ρₑ::Float64
    Ec::Float64
    Dₕₕ::Float64
    αₕₕ::Float64
    βₕₕ::Float64
    Cₕ::Float64
    δₕ::Float64
    r₀ₕₕ::Float64
    ϵₕ::Float64
    aₕ::Float64
    bₕ::Float64
    cₕ::Float64
    dₕ::Float64
    D_PdH::Float64
    α_PdH::Float64
    r₀PdH::Float64
    rₑ::Float64
    Ω::Float64
    fₑ::Float64

    atom_types::Vector{Symbol}
    cell::PeriodicCell
    cutoff::Float64
end

function PdH(atom_types::AbstractVector{Symbol}, cell::PeriodicCell, cutoff::Real)
    χ = 2.054085
    Φ = 0.216817
    δ = 8.414105
    β_PdPd = 7.221224
    η = 0.999999
    ρₑ = 3.316887

    # Iyad Hijazi, Yang Zhang & Robert Fuller (2018) A simple embeddedatom potential for Pd-H alloys, Molecular Simulation, 44:17, 1371-1379
    # Table 2
    Ec = 3.91 # Pd 

    # Iyad Hijazi, Yang Zhang & Robert Fuller (2018) A simple embeddedatom potential for Pd-H alloys, Molecular Simulation, 44:17, 1371-1379
    # Table 3
    Dₕₕ = 0.589510
    αₕₕ = 1.104827
    βₕₕ = 0.942490
    Cₕ = 2.145808
    δₕ = 0.942201
    r₀ₕₕ = 3.474173
    
    # Ciufo, R.A., Henkelman, G. Embedded atom method potential for hydrogen on palladium surfaces. J Mol Model 26, 336 (2020).
    ϵₕ = 0.054
    aₕ = 4.245745
    bₕ = 69.565833
    cₕ = 0.000104
    dₕ = 1.672605
    D_PdH = 0.123452
    α_PdH = 5.583967
    r₀PdH = 1.712318

    # Unit in angstrom and eV. assuming this is necessary
    rₑ = 2*1.37 # Metallic radius*2 https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
    Ω = 4/3π * rₑ^3 # Atomic volume
    fₑ = Ec / Ω

    PdH(χ, Φ, δ, β_PdPd, η, ρₑ, Ec, Dₕₕ, αₕₕ, βₕₕ, Cₕ, δₕ, r₀ₕₕ, ϵₕ, aₕ, bₕ, cₕ, dₕ, D_PdH, α_PdH, r₀PdH, rₑ, Ω, fₑ, atom_types, cell, cutoff)
end
        
function potential!(model::PdH, V::AbstractVector, R::AbstractMatrix)
    ########################
    # This has to be checked
    # They never say what it is
    ########################
    F(ρ) = -2.5
    F_Pd(ρ) = F(ρₑ) * (1 - model.η*log(ρ/model.ρₑ)) * (ρ/model.ρₑ)^model.η
    ρ_Pd(r) = model.fₑ * exp(-model.χ * (r - model.rₑ))

    ϕ_PdPd(r) = -model.Φ * (1 + model.δ*(r/model.rₑ - 1)) * exp(-model.β_PdPd*(r/model.rₑ-1))
    
    ϕₕₕ(r) = model.Dₕₕ * (model.βₕₕ * exp(-model.αₕₕ * (r - model.r₀ₕₕ)) - model.αₕₕ * exp(-model.βₕₕ * (r - model.r₀ₕₕ)))
    
    Fₕ(ρ) = -model.cₕ * (1/(2+model.dₕ)*(ρ+model.ϵₕ)^(2+model.dₕ) - (model.aₕ+model.bₕ)/(1+model.dₕ)*(ρ+model.ϵₕ)^(1+model.dₕ) + (model.aₕ*model.bₕ)/model.dₕ*(ρ+model.ϵₕ)^model.dₕ)
    ρₕ(r) = model.Cₕ * exp(-model.δₕ * r)

    ϕ_PdH(r) = model.D_PdH * (1 - exp(-model.α_PdH*(r - model.r₀PdH)))^2 - model.D_PdH
        
    function ϕ(r, i, j)
        if model.atom_types[i] == model.atom_types[j]
            if model.atom_types[i] == :Pd # Both Pd
                return ϕ_PdPd(r)
            else # Both H
                return ϕₕₕ(r)
            end
        else # One each
            return ϕ_PdH(r)
        end
    end
        
    pairs = get_pairs(R, model.cutoff, model.cell)
    E = 0.0
    densities = zeros(length(model.atom_types))
    for (i, j, R) in NeighbourLists.pairs(pairs)
        r = ustrip(u"Å", norm(R)u"bohr") # Convert distance to Å
        if model.atom_types[j] == :Pd
            densities[i] += ρ_Pd(r)
        else
            densities[i] += ρₕ(r)
        end
        E += ϕ(r, i, j)/2
    end
    
    for i=1:length(model.atom_types)
        if model.atom_types[i] == :Pd
            E += F_Pd(densities[i])
        else
            E += Fₕ(densities[i])
        end
    end
    V .= austrip(E*u"eV")
end
        
function derivative!(model::PdH, D::AbstractMatrix, R::AbstractMatrix)
    potential(x) = potential!(model, [0.0], x)[1]
    FiniteDiff.finite_difference_gradient!(D, potential, R)
end
