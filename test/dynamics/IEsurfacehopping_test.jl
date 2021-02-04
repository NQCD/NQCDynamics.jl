push!(LOAD_PATH, pwd())
using NonadiabaticMolecularDynamics
using NonadiabaticMolecularDynamics.InitialConditions
using LinearAlgebra
using Unitful
using UnitfulAtomic
using Revise
using Random, Distributions
using Plots

# IESH should proceed by the following steps:
#   1) initialize classical positions and momenta, initialize
#      electronic states, define surface on which to propagate
#   2) Do Verlet algorithm for the atomic positions, velocities electronic
#   3) Integrate the electronic wave function w/ 1e- Schroedinger equation
#   4) calculate switching probability, calculate if hop should occur
#   5) If no hop, GOTO 2; if hop, see if enough energy to do it, GOTO 2

n_states = 0# # number of electronic n_states
x = 5.0 # position of 1D particle
vinit = 0.005 # initial velocity
step = 0.01 # MD step
mass = 1.5 # atomic mass#
B_na = 1.0 # parameter that defines nonadiabatic coupling strength
n_DOF = 1
state = 1
ntrajes = 2

# 1a) Define adiabatic model, fom Model/scattering_anderson_holsteins.jl
#model = Models.ScatteringAndersonHolstein(N=n_states, a0=1, Bg=B_na, 
#Cg=1, alpha=-5.0, beta=.5)
#model = Models.TullyModelOne()
model = Models.Subotnik_A()
#R = rand(1,1)
#V = Hermitian(zeros(2,2))
#D = zero(R)
#D = [Hermitian(zeros(2, 2))]'
#Models.potential!(model,V,R)
#Models.derivative!(model,D,R)

# 1b) Initialize atomic parameters, i.e. moving atoms, of which there is one
atoms = Atoms{Float64}([:H])

# create Boltzmann distribution of momenta
temperature=300u"K"
#boltzmann = Normal(0.0, austrip(temperature)/atoms.masses[1])
boltzmann = Normal(0.0, austrip(temperature))

# create distribution of atomic positions from Monte Carlo sampling of the PESs
# for the number of trajectories and atoms
Δ = Dict([(:H, 0.3)])
monte_carlo = InitialConditions.MonteCarlo{Float64}(Δ, length(atoms), ntrajes, Int[])
sim = Simulation(atoms, model; DoFs=1, temperature=temperature)
outbolz = InitialConditions.run_monte_carlo_sampling(sim, monte_carlo, zeros(1, length(sim.atoms)))
#If sampled, gives an array of momenta and positions from the above initialized sampling
bolz_pos_momenta = PhasespaceDistribution(outbolz.R, boltzmann, (1,1))




###############################################################################
#
# Print the potential to look at adiabatic and diabatic states
#
##############################################################################
#points = [u for u in -10:0.1:30]
#p1 = zeros(size(points))
#p2 = zeros(size(points))
#p3 = zeros(size(points))
#eigs = zeros(length(points),model.n_states)
#diabats = zeros(length(points),model.n_states)
#for i=1:length(points)
#    p = zeros(1,1)
#    p[1] = points[i]
#    Models.potential!(model,V,p)
#    p1[i] = V[1,1]
#    p2[i] = V[2,2]
#    p3[i] = V[1,2]/0.0001
#    eigs[i,:] .=eigvals(V)
#end
# plot PESs
#plo=plot(points,p1, label="V1")
#plot!(points,p2, label="V2")
#plot!(points,eigs, label="eigs")
#plot!(points,diabats, label="diabats")

# plot hopping and coupling
#plo=plot(points,p3, label="Γ")
#plot!(points,eigs, label="eigs")
#xlabel!("x")
#ylabel!("Energy (a.u.)")
#savefig(plo,"Gammas.png")


##############################################################################
#
# Try dynamics
#
###############################################################################
#1c define the dynamics that will be used
# defined in main/Dynamics/iesh.jl
#IESH{T}(DoFs::Integer, atoms::Integer, states::Integer)
dynam = Dynamics.IESH{Float64}(1, 1, n_states+2)

# Initialize the simulation problem
# Simulation is defined in main/simulations.jl
#sim = Simulation(atoms, Models.TullyModelOne(), dynam; DoFs=1)
sim = Simulation(atoms, model, dynam; DoFs=1)
#calculate momentum

# Try using ensemble_problem
solution = Dynamics.run_ensemble(bolz_pos_momenta, (0.0,150000.0), sim; trajectories = ntrajes)
#println(solution.u[1])

#Do dynamics
#for i=1:ntrajes
#    fi = open("myfile_$i.txt", "w")
#    write(fi, "step, r (a.u.), energy (a.u.), p (a.u.), state\n")
#    p = fill(pall[i]*atoms.masses[1], sim.DoFs, length(sim.atoms))
#    r = fill(rall[i], sim.DoFs, length(sim.atoms)) 
#    println(p)

#    r = fill(-10., sim.DoFs, length(sim.atoms)) 
#    p = fill(vinit*atoms.masses[1], sim.DoFs, length(sim.atoms))
#    println(p)

    # intial state
    #k = rand(1:n_states+2)
#    k = 1
    #Initialize the surface hopping phasespace
    # postions, momenta, density matrix, state
    # ../../Dynamics/iesh.jl
#    z = SurfaceHoppingPhasespace(r,p, n_states+2, k)

    #display(k)

    #../../Dynamics/iesh.jl
    # Solution of Differentiall equations and propagation, step needs to be implemented
    #solution = Dynamics.run_trajectory(z, (0.0, 1000), sim, dt=step)
    # 
#    solution = Dynamics.run_trajectory(z, (0.0, 15000.0), sim)

    #    r1 = zeros(length(solution))
#    istate = zeros(length(solution))
#    r1 = [real(Dynamics.get_positions(u)[1,1]) for u in solution.u]
#    p2 = [real(Dynamics.get_momenta(u)[1,1]) for u in solution.u]
#    istate = [real(Dynamics.get_density_matrix(u)[1,1]) for u in solution.u]
#    eigs = zeros(length(solution),model.n_states)
#    p1 = zeros(size(solution))
#    V = Hermitian(zeros(2,2))
#    for j=1:length(solution)
#        p = zeros(1,1)
#        p[1] = r1[j]
#        Models.potential!(model,V,p)
#    p1[i] = V[1,1]
#    p2[i] = V[2,2]
#    p3[i] = V[1,2]/0.0001
#        eigs[j,:] .=eigvals(V)
#        p1[j] = eigs[j,1]
#        rrr = r1[j]
#        ppp = p1[j]
#        pppp = p2[j]
#        sss = istate[j]
#       write(fi, "$j $rrr $ppp $pppp $sss\n")
#    end
#    close(fi)
    
#println(solution.u[1])
    # for definitiion of u. and 
 #   a=plot([real(Dynamics.get_positions(u)[1,1]) for u in solution.u], [u.state for u in solution.u].-1, label="current surface", marker =2)
 #   plot!(r1, eigs[:,1]*100, label="Energy *100 (a.u.)", marker =2)
#plot(solution.t, [real(Dynamics.get_density_matrix(u)[1,1]) for u in solution.u], label="σ[1,1]", marker =2)
#plot!(solution.t, [real(Dynamics.get_density_matrix(u)[2,2]) for u in solution.u], label="σ[2,2]", marker =2)
#plot!(solution.t, [u.state for u in solution.u].-1, label="current surface")
#    xlabel!("x")
#    ylabel!("Energy (a.u.)")
#    display(plot(a))
#end
#savefig(plo,"traj_subA.png")

#display(plot(b))Plots