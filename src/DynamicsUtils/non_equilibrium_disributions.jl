using Unitful, UnitfulAtomic

# This script is called from `DynamicsUtils.jl` using the `include()` function



"""
    simple_non_eq(ϵ, μ, β, Δf_β1, Δf_β2)

Outputs a simple approximation to a non-equilibrium distribution, generated from finite temperature Fermi-Dirac distributions, that can be used for test purposes.
"""
function simple_non_eq(ϵ, μ, β, Δf_β1, Δf_β2)
    β != Inf || throw(error("Non-equilibrium distribution needs to be a finite temperature, β cannot be Inf."))
    Δf = DynamicsUtils.fermi.(ϵ, μ, Δf_β1) - fermi.(ϵ, μ, Δf_β2) #  create perterbation from 2 finite temperature Fermi-Dirac distributions
    F = DynamicsUtils.fermi.(ϵ, μ, β) + Δf  # add perterbation to a finite temperature Fermi-Dirac distribution
    return F
end




function integrate_ys(ys, bandmin, bandmax)
    h = (bandmax - bandmin)/length(ys)
    int = h*(ys[1] + ys[end])/2
    for Fi in ys[2:end-1]
        int = int + h*Fi
    end
    return int
end

function calc_dist_internal_energy(energies, distribution)
    bandwidth = energies[end] - energies[1]
    u = ((length(energies))/bandwidth) * integrate_ys(distribution.*energies, energies[1], energies[end])
    return u
end

function match_internal_energy(excitation_dist, energies, u_internal; c0=0.0, step= 0.0001, tol=0.001, maxcount=1e6, T_base=300)

    # # Internal energy of 0K fermi dirac distribution
    # T0 = 0.00001 # T = 0K approx
    # β0 = 1/austrip(T0*u"K")
    # FD0 = DynamicsUtils.fermi.(energies, 0.0, β0)
    # u0 = calc_dist_internal_energy(energies, FD0)
    
    # # Base distribution for neq dist
    T = T_base # K
    β = 1/austrip(T*u"K")
    global FD = DynamicsUtils.fermi.(energies, 0.0, β)

    dist_0 = FD .+ ( c0 .* (excitation_dist - reverse(excitation_dist)))
    global u = calc_dist_internal_energy(energies, dist_0)
    
    global c = copy(c0)
    global u_last = copy(u)
    global dist = copy(dist_0)

    count = 0
    # while (u - u0) != u_excitation
    while u  != u_internal
        count += 1

        # if isapprox((u - u0), u_internal, atol=u_internal*tol) == true
        if isapprox(u, u_internal, atol=u_internal*tol) == true
            # save previous iterations to the global variables
            global c = copy(c_i)
            global dist = copy(dist_i)
            global u_last = copy(u)
            
            break # break out of while loop
        elseif count >= maxcount
            throw(error("Max count reached, increase tolerance or vary step such that convergence occurs in given time."))
        else
            # change local variables for next trial
            global c_i = c0 + (count * step)
            global dist_i = FD .+ ( c_i .* (excitation_dist - reverse(excitation_dist)))
            global u = calc_dist_internal_energy(energies, dist_i)

        end
    end

    @info "Scaling factor optimised in $count iterations.\nOptimal value is $c."
    return dist
end



"""
    generate_NEQ_dist(excitation_energy::Float64, energy_grid, T_equivalent; T_ex = 300, T_base = 300)

Generates a non-equilibrium distribution, with an excitation up to a chosen energy represented as a step.
"""
function generate_NEQ_dist(excitation_energy::Float64, energy_grid, T_equivalent; T_ex=300, T_base=300)

    # step excitation
    exc_a = DynamicsUtils.fermi.(energy_grid, excitation_energy, 1/austrip(T_ex*u"K"))
    exc_b = reverse(DynamicsUtils.fermi.(energy_grid, 0.0, 1/austrip(T_base*u"K")))
    ex_dist = (exc_a + exc_b) .- 1 

    # calculate internal energy for the equivalent temperature Fermi Dirac distribution used in matching
    FD_equivalent = DynamicsUtils.fermi.(energy_grid, 0.0, 1/austrip(T_equivalent*u"K"))
    u_internal = calc_dist_internal_energy(energy_grid, FD_equivalent)

    # optimise NEQ excitaton amplitude to match require internal energy
    dist = match_internal_energy(ex_dist, energy_grid, u_internal; c0=0.0, step=0.00000001, tol=0.01, maxcount=1e9, T_base=T_base)

    return dist
end


