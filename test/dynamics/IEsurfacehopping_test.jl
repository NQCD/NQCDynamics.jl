push!(LOAD_PATH, pwd())
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using LinearAlgebra
using Unitful
using UnitfulAtomic
using Revise
using Random, Distributions
using Plots
using DelimitedFiles
using Profile
# using DataFrames

# IESH should proceed by the following steps:
#   1) initialize classical positions and momenta, initialize
#      electronic states, define surface on which to propagate
#   2) Do Verlet algorithm for the atomic positions, velocities electronic
#   3) Integrate the electronic wave function w/ 1e- Schroedinger equation
#   4) calculate switching probability, calculate if hop should occur
#   5) If no hop, GOTO 2; if hop, see if enough energy to do it, GOTO 2

ps = 98.22694788464062
Random.seed!(13)
#n_states = 202
n_states = 42
# # number of electronic n_states
vinit = 0.005 # initial velocity
mass = 1.5 # atomic mass#
B_na = 1.0 # parameter that defines nonadiabatic coupling strength
n_DOF = 1
ntrajes = 250
tspan = (0.0, 5.0u"ps")

# 1a) Define adiabatic model, fom Model/scattering_anderson_holsteins.jl
# model = Models.ScatteringAndersonHolstein(N=n_states, a0=1, Bg=B_na, 
# Cg=1, alpha=-5.0, beta=.5)
# model = Models.TullyModelOne()
# model = Models.Subotnik_A()
#model = Models.MiaoSubotnik(n_states=n_states, W=0.02, Gamma=0.0001)
model = Models.MiaoSubotnik(n_states=n_states, W=0.064, Gamma=0.0064)
# D = zero(R)
# D = [Hermitian(zeros(2, 2))]'
# Models.potential!(model,V,R)
# Models.derivative!(model,D,R)
# println(Models.potential!(model,V,R))
# Get eigenvalues
# Models.energy(model,R,n_states+2)

# 1b) Initialize atomic parameters, i.e. moving atoms, of which there is one
atoms = Atoms{Float64}([:H])



##############################################################################
#
# Try dynamics
#
###############################################################################

# create Boltzmann distribution of momenta
temperature = 300u"K"
temperature = 2.5*298u"K"
boltzmann = Normal(0.0, austrip(temperature) / atoms.masses[1])
boltzmann = Normal(0.0, austrip(temperature))

# create distribution of atomic positions from Monte Carlo sampling of the PESs
# for the number of trajectories and atoms
Δ = Dict([(:H, 0.3)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), ntrajes, Int[])
sim = Simulation(atoms, model; DoFs=1, temperature=temperature)
outbolz = InitialConditions.run_monte_carlo_sampling(sim, monte_carlo, zeros(1, length(sim.atoms)))
# If sampled, gives an array of momenta and positions from the above initialized sampling
bolz_pos_momenta = PhasespaceDistribution(outbolz.R, boltzmann, (1, 1))


# Save the energy
# Define the dynamics that will be used, see main/Dynamics/iesh.jl
dynam = Dynamics.IESH{Float64}(1, 1, n_states)
# dynam = Dynamics.FSSH{Float64}(1, 1, 2)


# Initialize the simulation problem; Simulation is defined in main/simulations.jl
# sim = Simulation(atoms, Models.TullyModelOne(), dynam; DoFs=1)
sim = Simulation(atoms, model, dynam; DoFs=1)


# Do dynamics
 for i=8:ntrajes
    nname=lpad(i,8,"0")
    r, p = rand(bolz_pos_momenta)
    p = -p * atoms.masses[1]
    println(i, r, p)

    # r = fill(-4.87, sim.DoFs, length(sim.atoms)) 
    #r = fill(0.13867834028319972, sim.DoFs, length(sim.atoms)) 
    # p = fill(vinit*atoms.masses[1], sim.DoFs, length(sim.atoms))
    #p = fill(2.7, sim.DoFs, length(sim.atoms))
    # println(r, p)

    # intial state (how to initialize a vector that shows which states are occupied)
    # For now, just initialize in ground state. Makes sense, anyhow, when the molecule
    # is hurling towards the surface, inititially
    # k = rand(1:n_states+2)

    k = round.(Int64, (zeros(n_states)))
    # k = zeros(n_states+2)
    l = Int(n_states / 2)
    k[1:l] .= 1
    println(k)
    # For FSSH
    # k = 1
    # Initialize the surface hopping phasespace
    # postions, momenta, density matrix, state-vector, see: ../../Dynamics/iesh.jl
    z = IESHPhasespace(r, p, n_states, k)
    


    # ../../Dynamics/callbacks.jl
    # and ../../Models/Models.jl
    # Solution of Differential equations and propagation, step needs to be implemented
    #cb, vals, vects = Dynamics.create_energy_saving_callback()
    #cb, vals = Dynamics.create_energy_saving_callback()
    cb, vals = Dynamics.create_impurity_saving_callback()
    # ../../Dynamics/iesh.jl
    @time solution = Dynamics.run_trajectory(z, (0.0, 100000.0), sim; callback=cb)
    # For testing only: w/o callback
    
    #@time solution = Dynamics.run_trajectory(z, (0.0, 12000000.0), sim)
    println("Finished")
    println(length(vals.t))

    open("trajectory_$nname.txt", "w") do fi
        write(fi, "step, r (a.u.), epot (a.u.), state, impurity population\n")
        outarray = zeros(length(vals.saveval), 5)
        for i = 1:length(vals.saveval)
            outarray[i,1] = vals.t[i]
            outarray[i,2] = vals.saveval[i][1]
            outarray[i,3] = vals.saveval[i][2]
            outarray[i,4] = vals.saveval[i][3]
            outarray[i,5] = vals.saveval[i][4]
        end
        writedlm(fi, outarray,'\t')
    end


# b = plot(outarray[:,1], outarray[:,2], label="energy", marker=2)
# xlabel!("x")
# ylabel!("Energy (a.u.)")

# a = plot(vals.t, outarray[:,4], label="Impurity pop", marker=2)
# xlabel!("time")
# ylabel!("Impurity Population (a.u.)")
# # # b = plot(solution.t[1:2:end], imp_pop, label="current surface", marker=2)
# # #    # plot!(outarray[:,2], outarray[:,3], label="current surface", marker =1)
# # # # a=plot(solution.t, [real(Dynamics.get_density_matrix(u)[20,20]) for u in solution.u], label="σ[1,1]", marker =2)
# # # # plot!(solution.t, [real(Dynamics.get_density_matrix(u)[21,21]) for u in solution.u], label="σ[2,2]", marker =2)
# # # # plot!(solution.t, [u.state[21] for u in solution.u].-1, label="current surface")
# display(plot(a))
# display(plot(b))

end
# # savefig(a,"traj_subA_G4Em4_W10G_test5.png")
# # savefig(b,"traj_imppop_G4Em4_W10G_test2.png")

# display(plot(b))