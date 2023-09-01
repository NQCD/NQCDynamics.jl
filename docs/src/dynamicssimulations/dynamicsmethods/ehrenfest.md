# [Ehrenfest molecular dynamics](@id ehrenfest-dynamics)

The Ehrenfest method is a mixed quantum-classical dynamics method in which the total wavefunction is factorized into slow (nuclear) variables, which are treated classically, and fast ones (electrons) which remain quantum-mechanical. In the Ehrenfest method, nuclei move according to classical mechanics on a potential energy surface given by the expectation value of the electronic Hamiltonian. 

The time dependence of the electronic wavefunction is expanded into an adiabatic basis and follows the time-dependent Schr\"odinger equation.
```math
i\hbar \dot{c}_i(t) = V_i(\mathbf{R}) c_i (t)
- i\hbar \sum_j \dot{\mathbf{R}} \cdot \mathbf{d}_{ij}(\mathbf{R})c_j(t)
```

## Example
Below the example of the Ehrenfest implementation is presented, using model from [Ananth2007](@cite).

At the start, we assign `atoms` and initialise the simulation using the mass and model from [NQCModels.jl](@ref).
```@example ehrenfest
using NQCDynamics

atoms = Atoms(1980)
sim = Simulation{Ehrenfest}(atoms, AnanthModelOne())
```

Next, the initial distribution is created:
```@example ehrenfest
using Distributions
e = 0.03
k = sqrt(e*2*atoms.masses[1])
r = Normal(-5, 1/sqrt(0.25))
v = k / atoms.masses[1]
distribution = DynamicalDistribution(v, r, size(sim))* PureState(1, Adiabatic())
```

To run an ensemble simulation we additionally choose number of trajectories `n_traj` and timespan `tspan` and we pass all the established settings to the `run_dynamics` function.
In this example we output velocities by specifying `output=OutputVelocity` and store the final values in the `final_velocities` array. Following that, we calculate final momenta.
```@example ehrenfest
n_traj = 10
tspan = (0.0, 3000.0)
solution = run_dynamics(sim, tspan, distribution; 
    trajectories=n_traj, output=OutputVelocity, dt=1.0)
final_velocities = [r[:OutputVelocity][end] for r in solution]
momenta = reduce(vcat, final_velocities*atoms.masses[1])
```

```@example ehrenfest
using Plots
histogram(momenta)
xlims!(-20,20)
```
