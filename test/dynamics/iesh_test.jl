push!(LOAD_PATH, pwd())
using Test
using NonadiabaticMolecularDynamics
using OrdinaryDiffEq

# Test setup of phasespace?
@test Dynamics.IESH{Float64}(3,3,12) isa Dynamics.IESH
atoms = Atoms{Float64}([:H])
n_states = 12 # number of states in the system
model = Models.MiaoSubotnik(n_states=n_states, W=0.064, Gamma=0.0064)
dynam = Dynamics.IESH{Float64}(1, 1, n_states)
sim = Simulation(atoms, model, dynam; DoFs=1)

R = zeros(sim.DoFs, length(sim.atoms))
P = rand(sim.DoFs, length(sim.atoms))
P = fill(2.7, sim.DoFs, length(sim.atoms))
#println(P)
#P = [0.8]

# Vector of initial populations
k = round.(Int64, (zeros(n_states)))
l = Int(n_states / 2)
k[1:l] .= 1
u = IESHPhasespace(R, P, n_states, k)
du = zero(u)
# set some velocity
du.x.x[1] .= P

# # test if states on diagonal elements
# @test tr(Dynamics.get_density_matrix(u)) ≈ Complex.(l)

# test if complex, diagonal matrix
A = zeros(Complex,length(k), length(k))
for i=1:l
    A[i,i] = 1
end
@test Dynamics.get_density_matrix(u) ≈ A

#Test if nonadiabatic coupling behave reasonably
Dynamics.motion!(du, u, sim, 1.0)
Dynamics.evaluate_nonadiabatic_coupling!(sim)
@test sim.method.nonadiabatic_coupling ≈ -sim.method.nonadiabatic_coupling'
##

# Test new state generatation
problem = ODEProblem(Dynamics.motion!, u, (0.0, 1.0), sim)
integrator = init(problem, Tsit5())
# necessary to modifiy dt outside of init().
set_proposed_dt!(integrator, 1.0)
# Set high values for density matrix and nonadiabatic coupling to force hopping
# Hopping is forced from the highest occupied to the lowest unoccupied
integrator.u.x.x[3][l+1,l] = 100
sim.method.nonadiabatic_coupling[1][l,l+1] = 100
Dynamics.update_hopping_probability!(integrator)

# Test that new states are correctly selected and hopping prob.
# is above 0
@testset "New states" begin
    @test 0 < integrator.p.method.hopping_probability[1]
    @test l == integrator.p.method.hopping_probability[2]
    @test l + 1 == integrator.p.method.hopping_probability[3]
end
##

# Check momentum rescaling
@testset "Hopping" begin
    # Check that hop doesn't happened (small momentum)
    # Read out state vector
    k = integrator.u.state
    ktest = copy(k)
    Dynamics.affect!(integrator)
    @test ktest == integrator.u.state

    # Check that hop happens with large moment
    # Set larger velocity
    get_positions(get_du(integrator)) .= 10.
    # Build expected state vector
    ktest[l] = 0
    ktest[l+1] = 1
    Dynamics.affect!(integrator)
    @test ktest == integrator.u.state
end
##

# Check momentum rescaling happening with low momentum
#get_momenta(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
#Dynamics.motion!(get_du(integrator), integrator.u, sim, 0.0) # Update velocity
#cont = Dynamics.calculate_rescaling_constant!(integrator, k)
#println(cont)
#@test cont == false # Check hop is rejected 




# set such that hopping can't happen
