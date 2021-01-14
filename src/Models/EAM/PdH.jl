export PdH

using LinearAlgebra
using Unitful
using FiniteDiff
using ForwardDiff

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

########################
# This has to be checked
# They never say what it is
########################
F(ρ) = -2.5
F_Pd(ρ, ρₑ, η) = F(ρₑ) * (1 - η*log(ρ/ρₑ)) * (ρ/ρₑ)^η

ρ_Pd(r, fₑ, χ, rₑ) = fₑ * exp(-χ * (r - rₑ))

ϕ_PdPd(r, Φ, δ, rₑ, β_PdPd) = -Φ * (1 + δ*(r/rₑ - 1)) * exp(-β_PdPd*(r/rₑ-1))

ϕₕₕ(r, Dₕₕ, βₕₕ, αₕₕ, r₀ₕₕ) = Dₕₕ * (βₕₕ * exp(-αₕₕ * (r - r₀ₕₕ)) - αₕₕ * exp(-βₕₕ * (r - r₀ₕₕ)))

Fₕ(ρ, aₕ, bₕ, cₕ, dₕ, ϵₕ) = -cₕ * (1/(2+dₕ)*(ρ+ϵₕ)^(2+dₕ) - (aₕ+bₕ)/(1+dₕ)*(ρ+ϵₕ)^(1+dₕ) + (aₕ*bₕ)/dₕ*(ρ+ϵₕ)^dₕ)

ρₕ(r, Cₕ, δₕ) = Cₕ * exp(-δₕ * r)

ϕ_PdH(r, D_PdH, α_PdH, r₀PdH) = D_PdH * (1 - exp(-α_PdH*(r - r₀PdH)))^2 - D_PdH
        
function potential!(model::PdH, V::AbstractVector, R::AbstractMatrix)
        
    function ϕ(r, i, j)
        if model.atom_types[i] == model.atom_types[j]
            if model.atom_types[i] == :Pd # Both Pd
                return ϕ_PdPd(r, model.Φ, model.δ, model.rₑ, model.β_PdPd)
            else # Both H
                return ϕₕₕ(r, model.Dₕₕ, model.βₕₕ, model.αₕₕ, model.r₀ₕₕ)
            end
        else # One each
            return ϕ_PdH(r, model.D_PdH, model.α_PdH, model.r₀PdH)
        end
    end
        
    pairs = get_pairs(R, model.cutoff, model.cell)
    E = 0.0
    densities = zeros(eltype(R), length(model.atom_types))
    for (i, j, R) in NeighbourLists.pairs(pairs)
        r = ustrip(auconvert(u"Å", norm(R))) # Convert distance to Å
        if model.atom_types[j] == :Pd
            densities[i] += ρ_Pd(r, model.fₑ, model.χ, model.rₑ)
        else
            densities[i] += ρₕ(r, model.Cₕ, model.δₕ)
        end
        E += ϕ(r, i, j)/2
    end
    
    for i=1:length(model.atom_types)
        if model.atom_types[i] == :Pd
            E += F_Pd(densities[i], model.ρₑ, model.η)
        else
            E += Fₕ(densities[i], model.aₕ, model.bₕ, model.cₕ, model.dₕ, model.ϵₕ)
        end
    end
    V .= austrip(E*u"eV")
end

function derivative!(model::PdH, D::AbstractMatrix, R::AbstractMatrix)
    potential(x) = potential!(model, zeros(eltype(x), 1), x)[1]
    D .= ForwardDiff.gradient(potential, R)
end
