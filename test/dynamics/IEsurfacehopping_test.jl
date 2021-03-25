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
n_states = 6
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
# model = Models.TullyModelOne()
# model = Models.Subotnik_A()
model = Models.MiaoSubotnik(n_states=n_states, W=0.004, Gamma=0.0004)
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
# for i=1:ntrajes
#    nname=lpad(i,8,"0")
#    fi = open("trajectory_$nname.txt", "w")
#    write(fi, "step, r (a.u.), p (a.u.),ekin (a.u.), energy (a.u.), state\n")
r, p = rand(bolz_pos_momenta)
p = -p * atoms.masses[1]
println(r, p)

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
    # z = SurfaceHoppingPhasespace(r,p, 2, k)


    # display(k)

    # ../../Dynamics/callbacks.jl
    # and ../../Models/Models.jl
    # Solution of Differential equations and propagation, step needs to be implemented
    # HUOM!!!! create_energy_saving_callback pitää korjata! Energy on väärin.

    # This will get more complicated, since the structure will need to change.
    # What I need are both the eigenvalues and the eigenvectors (which could, 
    #in principle be encoded in an n_state+1 nstate matrix where the last line 
    # contains the eigenvalues)
    # It's a wonder it works at all, that with its reference to SurfaceHoppingPhasespace.
    #cb, vals, vects = Dynamics.create_energy_saving_callback()
    #cb, vals = Dynamics.create_energy_saving_callback()
    cb, vals = Dynamics.create_eigen_saving_callback()
    # ../../Dynamics/iesh.jl
@time solution = Dynamics.run_trajectory(z, (0.0, 10.0), sim; callback=cb)
    # For testing only: w/o callback
    
    # solution = Dynamics.run_trajectory(z, (0.0, 15000.0), sim)
println("Finished")

println(size(solution.t))
println(vals.saveval)
    
    # array for output
# outarray = zeros(length(solution.t[1:2:end]), 6)
# outarray[:,1] = solution.t[1:2:end]
# outarray[:,2] = [real(Dynamics.get_positions(u)[1,1]) for u in solution.u][1:2:end]
# outarray[:,3] = [real(Dynamics.get_momenta(u)[1,1]) for u in solution.u][1:2:end]
# outarray[:,6] = [u.state[3] for u in solution.u][1:2:end] .- 1
# outarray[:,4] = NonadiabaticMolecularDynamics.evaluate_kinetic_energy.(Ref(sim), get_momenta.(solution.u))[1:2:end]
# outarray[:,5] = vals.saveval

# #Calculate the impurity population by turning the density matrix into its diabatic
# # represenation using the eigenvectors.
# imp_pop = zeros(length(solution.u))
# g = Int(n_states/2)
# mat_j = zeros(g, g)
# mat_ρ = zeros(ComplexF64,g, g)
# i = 0
# s = 0 
# for u in solution.u
#     #println(det(u.x.x[3]))
#     #println(real(det(u.x.x[3]))^2+ imag(det(u.x.x[3]))^2 )
#     i = 0
#     global s = s + 1
#     for n=1:n_states
#         j = 0
#         if (u.state[n] ==1)
#             i = i + 1
#             mat_j[i,i] = u.state[n]
#             for m = 1: n_states
#                 if (u.state[m] ==1)
#                     j = j +1
#                     mat_ρ[i,j] = u.x.x[3][n,m]
#                 end
#             end 
#         end
#         imp_pop[s] = imp_pop[s] + (1-real(u.x.x[3][n,n]))^2*u.state[n]
#     end
#     println(real(det(mat_j*mat_ρ))^2+imag(det(mat_j*mat_ρ))^2)
#     #imp_pop[i] = 1 - real(det(u.x.x[3]))^2+ imag(det(u.x.x[3]))^2  #imp_pop[i]
#     #imp_pop[i] = imp_pop[i]
#     #imp_pop[s] = real(det(mat_j*mat_ρ))^2+imag(det(mat_j*mat_ρ))^2
# end 
    

#         # Needs a file to write things to
#     # for j=1:length(solution.t[1:2:end])
#     #    writedlm(fi,[ outarray[j,1] outarray[j,2] outarray[j,3] outarray[j,4] outarray[j,5] outarray[j,6]])
#     # end
    
#    # println(outarray[:,3])
#    # println(outarray[:,2])
#    # a=plot([real(Dynamics.get_positions(u)[1,1]) for u in solution.u], [u.state[20] for u in solution.u].-1, label="current surface", marker =2)
# a = plot(outarray[:,2], outarray[:,5] * 100, label="current surface", marker=2)
# b = plot(solution.t, imp_pop, label="current surface", marker=2)
#    # plot!(outarray[:,2], outarray[:,3], label="current surface", marker =1)
# # a=plot(solution.t, [real(Dynamics.get_density_matrix(u)[20,20]) for u in solution.u], label="σ[1,1]", marker =2)
# # plot!(solution.t, [real(Dynamics.get_density_matrix(u)[21,21]) for u in solution.u], label="σ[2,2]", marker =2)
# # plot!(solution.t, [u.state[21] for u in solution.u].-1, label="current surface")
# xlabel!("x")
# ylabel!("Energy (a.u.)")
# display(plot(a))
# display(plot(b))

# # end
# savefig(a,"traj_subA_G4Em4_W10G_test5.png")
# savefig(b,"traj_imppop_G4Em4_W10G_test2.png")

# display(plot(b))Plots