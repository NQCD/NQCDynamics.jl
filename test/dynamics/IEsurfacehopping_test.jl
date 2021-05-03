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
using BenchmarkTools


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
ntrajes = 1
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
# temperature = 300u"K"
# temperature = 2.5*298u"K"
# boltzmann = Normal(0.0, austrip(temperature) / atoms.masses[1])
# boltzmann = Normal(0.0, austrip(temperature))

# create distribution of atomic positions from Monte Carlo sampling of the PESs
# for the number of trajectories and atoms
### SOMETHING HERE IS BROKEN
# Δ = Dict([(:H, 0.3)])
# tsim = Simulation(atoms, model; DoFs=1, temperature=temperature)
# outbolz = InitialConditions.run_monte_carlo_sampling(tsim, zeros(1, length(tsim.atoms)), Δ, ntrajes)
# # If sampled, gives an array of momenta and positions from the above initialized sampling
# bolz_pos_momenta = DynamicalDistribution(boltzmann, outbolz.R, (1, 1))




# # Initialize the simulation problem; Simulation is defined in src/simulation_constructors.jl 
sim = Simulation{IESH}(atoms, model; DoFs=1)


# # Do dynamics
  for i=1:ntrajes
     nname=lpad(i,8,"0")
#     r, v = rand(bolz_pos_momenta)
#     v = -v
#     println(i," ", r, " ", v)

r = fill(-5.0, sim.DoFs, length(sim.atoms))
v = fill(8.9, sim.DoFs, length(sim.atoms)) ./ sim.atoms.masses[1]
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
    z = SurfaceHoppingVariables(v, r, n_states, k)
    


    #run_trajectory defined in src/Dynamics.jl
    # callbacks defined in: src/Dynamics/callbacks.jl and src/Models/Models.jl
#    @time solution = Dynamics.run_trajectory(z, (0.0, 100.0), sim; output=(:density_matrix, :state))
    # Save impurity does not seem to be defined at the moment?
    #@time solution = Dynamics.run_trajectory(z, (0.0, 5000.0), sim; output=(:save_impurity))
    @time solution = Dynamics.run_trajectory(z, (0.0, 100000.0), sim, dt=1, 
                                             adaptive=false; output=(:position, :save_impurity))
    
    #@time solution = Dynamics.run_trajectory(z, (0.0, 12000000.0), sim)
    println("Finished")
    println(length(solution.t))

    #open("trajectory_$nname.txt", "w") do fi
    #    write(fi, "step, r (a.u.), epot (a.u.), state, impurity population\n")
        outarray = zeros(length(solution), 5)
        for i = 1:length(solution)
            outarray[i,1] = solution.t[i]
            outarray[i,2] = solution.position[i][1]
            outarray[i,3] = solution.save_impurity[i][2]
            # outarray[i,2] = solution.save_impurity[i][1]
            # outarray[i,3] = solution.save_impurity[i][2]
            # outarray[i,4] = solution.save_impurity[i][3]
            # outarray    [i,5] = solution.save_impurity[i][4]
        end
    #    writedlm(fi, outarray,'\t')
    #end


b = plot(outarray[:,2], outarray[:,3], label="energy", marker=2)
xlabel!("x")
ylabel!("Energy (a.u.)")

# # # a = plot(vals.t, outarray[:,4], label="Impurity pop", marker=2)
# # # xlabel!("time")
# # # ylabel!("Impurity Population (a.u.)")
# # # # # b = plot(solution.t[1:2:end], imp_pop, label="current surface", marker=2)
# # # # #    # plot!(outarray[:,2], outarray[:,3], label="current surface", marker =1)
# # # # # # a=plot(solution.t, [real(Dynamics.get_density_matrix(u)[20,20]) for u in solution.u], label="σ[1,1]", marker =2)
# # # # # # plot!(solution.t, [real(Dynamics.get_density_matrix(u)[21,21]) for u in solution.u], label="σ[2,2]", marker =2)
# # # # # # plot!(solution.t, [u.state[21] for u in solution.u].-1, label="current surface")
# # # display(plot(a))
display(plot(b))

end
# # # savefig(a,"traj_subA_G4Em4_W10G_test5.png")
# # # savefig(b,"traj_imppop_G4Em4_W10G_test2.png")

#display(plot(b))