export PdH

using LinearAlgebra
using Unitful
using FiniteDifferences

"""
Ciufo, R.A., Henkelman, G. Embedded atom method potential for hydrogen on palladium surfaces. J Mol Model 26, 336 (2020).
"""
struct PdH <: Model

    n_states::UInt
    get_V0::Function
    get_D0::Function

    # a0 = 3.89 Å
    function PdH(atom_types::Vector{Symbol}, cell::Atoms.AbstractCell, cutoff::AbstractFloat)

        # Iyad Hijazi, Yang Zhang & Robert Fuller (2018) A simple embeddedatom potential for Pd-H alloys, Molecular Simulation, 44:17, 1371-1379
        # Table 1
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
        
        ########################
        # This has to be checked
        # They never say what it is
        ########################
        F(ρ) = -2.5
        F_Pd(ρ) = F(ρₑ) * (1 - η*log(ρ/ρₑ)) * (ρ/ρₑ)^η
        fₑ = Ec / Ω
        ρ_Pd(r) = fₑ * exp(-χ * (r - rₑ))

        ϕ_PdPd(r) = -Φ * (1 + δ*(r/rₑ - 1)) * exp(-β_PdPd*(r/rₑ-1))
        
        ϕₕₕ(r) = Dₕₕ * (βₕₕ * exp(-αₕₕ * (r - r₀ₕₕ)) - αₕₕ * exp(βₕₕ * (r - r₀ₕₕ)))
        
        Fₕ(ρ) = -cₕ * (1/(2+dₕ)*(ρ+ϵₕ)^(2+dₕ) - (aₕ+bₕ)/(1+dₕ)*(ρ+ϵₕ)^(1+dₕ) + (aₕ*bₕ)/dₕ*(ρ+ϵₕ)^dₕ)
        ρₕ(r) = Cₕ * exp(-δₕ * r)

        ϕ_PdH(r) = D_PdH * (1 - exp(-α_PdH*(r - r₀PdH)))^2 - D_PdH
        
        sum_all_except(a, i) = sum(a[1:end .!= i])
        
        function ϕ(r, i, j)
            if atom_types[i] == atom_types[j]
                if atom_types[i] == :Pd # Both Pd
                    return ϕ_PdPd(r)
                else # Both H
                    return ϕₕₕ(r)
                end
            else # One each
                return ϕ_PdH(r)
            end
        end
        
        function get_V0(R)
            pairs = get_pairs(R, cutoff, cell)
            E = 0.0
            densities = zeros(length(atom_types))
            for (i, j, R) in NeighbourLists.pairs(pairs)
                r = ustrip(u"Å", norm(R)u"bohr") # Convert distance to Å
                if atom_types[j] == :Pd
                    densities[i] += ρ_Pd(r)
                else
                    densities[i] += ρₕ(r)
                end
                E += ϕ(r, i, j)/2
            end
            
            for i=1:length(atom_types)
                if atom_types[i] == :Pd
                    E += F_Pd(densities[i])
                else
                    E += Fₕ(densities[i])
                end
            end
            austrip(E*u"eV")
        end
        
        # get_D0(R) = ForwardDiff.gradient(get_V0, R) # ForwardDiff autodiff, error with NeighbourLists
        # get_D0(R) = get_V0'(R) # Zygote autodiff, crashes
        get_D0(R) = grad(central_fdm(5, 1), get_V0, R)[1] # First derivative, 5th order

        new(1, get_V0, get_D0)
    end
end
