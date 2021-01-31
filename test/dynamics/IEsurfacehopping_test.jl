push!(LOAD_PATH, pwd())
using NonadiabaticMolecularDynamics
using LinearAlgebra: Hermitian
using Unitful
using Revise
using Random
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
vinit = -0.005 # initial velocity
step = 0.01 # MD step
mass = 1.5 # atomic mass#
B_na = 1.0 # parameter that defines nonadiabatic coupling strength
n_DOF = 1
state = 1


# 1a) Define adiabatic model, fom Model/scattering_anderson_holsteins.jl
#model = Models.ScatteringAndersonHolstein(N=n_states, a0=1, Bg=B_na, 
#Cg=1, alpha=-5.0, beta=.5)
#model = Models.TullyModelOne()
model = Models.Subotnik_A()
R = rand(1,1)
V = Hermitian(zeros(2,2))
Models.potential!(model,V,R)
println(V)

points = [u for u in -10:0.1:30]

p1 = zeros(size(points))
p2 = zeros(size(points))

for i=1:length(points)
    R = zeros(1,1)
    R[1] = points[i]
    Models.potential!(model,V,R)
    p1[i] = V[1,1]
    p2[i] = V[2,2]
    #println(V)

end

plot(points,model)

# 1b) Initialize atomic parameters, i.e. moving atoms, of which there is one
#atoms = Atoms{Float64}([:H])

#1c define the dynamics that will be used
# defined in main/Dynamics/iesh.jl
#IESH{T}(DoFs::Integer, atoms::Integer, states::Integer)
#dynam = Dynamics.IESH{Float64}(1, 1, n_states+2)

# Initialize the simulation problem
# Simulation is defined in main/simulations.jl
#sim = Simulation(atoms, Models.TullyModelOne(), dynam; DoFs=1)
#sim = Simulation(atoms, model, dynam; DoFs=1)
#calculate momentum
#r = fill(x, sim.DoFs, length(sim.atoms)) 
#p = fill(vinit*atoms.masses[1], sim.DoFs, length(sim.atoms))

# intial state
#k = rand(1:n_states+2)
#k = 1
#Initialize the surface hopping phasespace
# postions, momenta, density matrix, state
# ../../Dynamics/iesh.jl
#z = SurfaceHoppingPhasespace(r,p, n_states+2, k)
#display(k)

#../../Dynamics/iesh.jl
# Solution of Differentiall equations and propagation
#solution = Dynamics.run_trajectory(z, (0.0, 1000), sim, dt=step)
# Problem is presumably in updat_electronics
#solution = Dynamics.run_trajectory(z, (0.0, 2500.0), sim)


#a =plot(solution.t, [real(Dynamics.get_positions(u)[1,1]) for u in solution.u], label="r")
#plot(solution.t, [real(Dynamics.get_density_matrix(u)[1,1]) for u in solution.u], label="σ[1,1]", marker =2)
#plot!(solution.t, [real(Dynamics.get_density_matrix(u)[2,2]) for u in solution.u], label="σ[2,2]", marker =2)
#plot!(solution.t, [u.state for u in solution.u].-1, label="current surface")
#display(plot(a))
#display(plot(b))