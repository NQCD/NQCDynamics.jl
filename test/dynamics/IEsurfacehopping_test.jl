push!(LOAD_PATH, pwd())
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using LinearAlgebra
using Unitful
using UnitfulAtomic
using Revise
using Random, Distributions
using Plots
#using DelimitedFiles
#using Profile
# using DataFrames
#using BenchmarkTools


# IESH should proceed by the following steps:
#   1) initialize classical positions and momenta, initialize
#      electronic states, define surface on which to propagate
#   2) Do Verlet algorithm for the atomic positions, velocities electronic
#   3) Integrate the electronic wave function w/ 1e- Schroedinger equationa
#   4) calculate switching probability, calculate if hop should occur
#   5) If no hop, GOTO 2; if hop, see if enough energy to do it, GOTO 2

ps = 98.22694788464062
Random.seed!(13)
#n_states = 202
n_states = 4
#n_states = 502
# # number of electronic n_states
vinit = 0.005 # initial velocity
mass = 1.5 # atomic mass#
B_na = 1.0 # parameter that defines nonadiabatic coupling strength
n_DOF = 1
#ntrajes = 1
tspan = (0.0, 5.0u"ps")

# 1a) Define adiabatic model, fom Model/scattering_anderson_holsteins.jl
# model = Models.ScatteringAndersonHolstein(N=n_states, a0=1, Bg=B_na, 
# Cg=1, alpha=-5.0, beta=.5)
#model = Models.MiaoSubotnik(n_states=n_states, W=0.064, Gamma=0.0064)
model = Models.MiaoSubotnik(n_states=n_states, W=0.0064, Gamma=0.064)
#model = Models.TullyModelOne()
#println(model.n_states)

# 1b) Initialize atomic parameters, i.e. moving atoms, of which there is one
atoms = Atoms{Float64}([:H])



##############################################################################
#
# Try dynamics
#
###############################################################################

# create Boltzmann distribution of momenta
temperature = 2.5*298u"K"

# create distribution of atomic positions from Monte Carlo sampling of the PESs
# for the number of trajectories and atoms
# atoms and step size
Δ = Dict([(:H, 0.3)])
tsim = Simulation{Classical}(atoms, model; DoFs=1, temperature=temperature)
R0 = zeros(3, length(tsim.atoms))
Δ = Dict([(:H, 0.3)])
R0 = zeros(1, length(tsim.atoms))
passes = 1000
outbolz = MetropolisHastings.run_monte_carlo_sampling(tsim, R0, Δ, passes)


# If sampled, gives an array of momenta and positions from the above initialized sampling
#bolz_pos_momenta = DynamicalDistribution(boltzmann, outbolz.R, (1, 1))




# # Initialize the simulation problem; Simulation is defined in src/simulation_constructors.jl 
sim = Simulation{wave_IESH}(atoms, model; DoFs=1)


# # Do dynamics
#for i=1:ntrajes
i = 1
     nname=lpad(i,8,"0")
#     r, v = rand(bolz_pos_momenta)
#     v = -v
#     println(i," ", r, " ", v)

    #r = fill(-5.0, sim.DoFs, length(sim.atoms))
    r = rand(outbolz.R)
    v = fill(sqrt(rand(outbolz.energy)*20/sim.atoms.masses[1]), sim.DoFs, length(sim.atoms))
    #v = fill(5/sim.atoms.masses[1], sim.DoFs, length(sim.atoms))
    println(r, v)

    # intial state (how to initialize a vector that shows which states are occupied)
    # For now, just initialize in ground state. Makes sense, anyhow, when the molecule
    # is hurling towards the surface, inititially
    # This is defined in src/Dynamics/SurfaceHopping/surface_hopping_variables.jl
    # as: SurfaceHoppingVariables(v::AbstractArray, r::AbstractArray, n_states::Integer, state::Integer)
    k = round.(Int64, (zeros(n_states)))
    l = Int(n_states / 2)
    k[1:l] .= 1
    println(k)
    z = SurfaceHoppingVariablesIESH(v, r, n_states, k)
    


    #run_trajectory defined in src/Dynamics.jl
    # callbacks defined in: src/Dynamics/callbacks.jl and src/Models/Models.jl
    # Save impurity does not seem to be defined at the moment?
    #@time solution = Dynamics.run_trajectory(z, (0.0, 10000000.0), sim; output=(:save_impurity))
    @time solution = Dynamics.run_trajectory(z, (0.0, 20000.0), sim, dt=10, adaptive=false; 
                                             output=(:position, :save_impurity))
    println("Finished")
    println(length(solution.t))

    #open("trajectory_$nname.txt", "w") do fi
    #    write(fi, "step, r (a.u.), epot (a.u.), state, impurity population\n")
    aka = Int(ceil(length(solution)/10))-1
    #aka = Int(length(solution))
    j = 0
    println(aka)
        outarray = zeros(aka, 5)
        for i = 1:length(solution)
            if (mod(i,10)==0)
                global j = j + 1
                outarray[j,1] = solution.t[i]
                outarray[j,2] = solution.position[i][1]
                outarray[j,3] = solution.save_impurity[i][2]
                outarray[j,4] = solution.save_impurity[i][3]
                outarray[j,5] = solution.save_impurity[i][4]
            # outarray[i,3] = solution.save_impurity[i][2]
            # outarray[i,4] = solution.save_impurity[i][3]
            # outarray    [i,5] = solution.save_impurity[i][4]
            end
        end
    #    writedlm(fi, outarray,'\t')
    #end


b = plot(outarray[:,1], outarray[:,2], label="energy", marker=2)
xlabel!("t")
ylabel!("x")
a = plot(outarray[:,2], outarray[:,3], label="energy", marker=2)
 xlabel!("x")
 ylabel!("Energy (a.u.)")
# a = plot(outarray[:,1], outarray[:,5], label="energy", marker=2)
# xlabel!("time")
# ylabel!("Impurity Population")

display(plot(b))
display(plot(a))

#end
# # # savefig(a,"traj_subA_G4Em4_W10G_test5.png")
# # # savefig(b,"traj_imppop_G4Em4_W10G_test2.png")

#display(plot(b))