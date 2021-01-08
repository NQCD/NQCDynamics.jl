using Test
using NonadiabaticMolecularDynamics
using DifferentialEquations

@test Dynamics.FSSH{Float64}(3, 3, 2) isa Dynamics.FSSH
atoms = Atoms{Float64}([:H])
sim = Simulation(atoms, Models.DoubleWell(), Dynamics.FSSH{Float64}(1, 1, 2); DoFs=1)

R = zeros(sim.DoFs, length(sim.atoms)) 
P = rand(sim.DoFs, length(sim.atoms)) 
u = SurfaceHoppingPhasespace(R, P, 2, 1)
du = zero(u)

@test Dynamics.get_density_matrix(u) ≈ Complex.([1 0; 0 0])

Dynamics.motion!(du, u, sim, 0.0)

Dynamics.evaluate_nonadiabatic_coupling!(sim)
@test sim.method.nonadiabatic_coupling ≈ -sim.method.nonadiabatic_coupling'

problem = ODEProblem(Dynamics.motion!, u, (0.0, 1.0), sim)
integrator = init(problem, Tsit5(), callback=Dynamics.fssh_callback)
Dynamics.update_hopping_probability!(integrator)

@testset "select_new_state" begin
    @test 3 == Dynamics.select_new_state([0, 0, 1.0], 2)
    @test 0 == Dynamics.select_new_state([0, 1.0, 0.0], 2)
    @test 1 == Dynamics.select_new_state([1.0, 0.0, 0.0], 2)
end

@testset "calculate_rescaling_constant" begin
    get_momenta(integrator.u) .= 0.0 # Set momentum to zero to force frustrated hop
    Dynamics.motion!(get_du(integrator), integrator.u, sim, 0.0) # Update velocity
    cont = Dynamics.calculate_rescaling_constant!(integrator, 2)
    @test cont == false # Check hop is rejected 

    get_momenta(integrator.u) .= 1e5 # Set momentum to big to force hop
    Dynamics.motion!(get_du(integrator), integrator.u, sim, 0.0) # Update velocity
    cont = Dynamics.calculate_rescaling_constant!(integrator, 2)
    @test cont == true # Check hop is accepted

    @test 1.2 ≈ Dynamics.calculate_potential_energy_change([0.0, 0.1, 0.9, -0.3], 3, 4)
end

@testset "execute_hop!" begin
    get_momenta(integrator.u) .= 1e4 # Set high momentum to ensure successful hop
    Dynamics.motion!(get_du(integrator), integrator.u, sim, 0.0)
    ΔE = Dynamics.calculate_potential_energy_change(integrator.p.calculator.eigenvalues, 2, 1)

    cont = Dynamics.calculate_rescaling_constant!(integrator, 2)
    @test cont == true

    KE_initial = sum(get_momenta(integrator.u).^2 ./ integrator.p.atoms.masses')/2
    @test integrator.u.state == 1
    Dynamics.execute_hop!(integrator, 2)
    @test integrator.u.state == 2 # Check state has changed
    KE_final = sum(get_momenta(integrator.u).^2 ./ integrator.p.atoms.masses')/2
    @test KE_final - KE_initial ≈ -ΔE # Test for energy conservation
end
