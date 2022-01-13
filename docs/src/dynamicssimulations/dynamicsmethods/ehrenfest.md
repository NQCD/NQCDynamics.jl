# [Ehrenfest molecular dynamics](@id ehrenfest-dynamics)

The Ehrenfest method is a mixed quantum-classical dynamics (MQCD) method in which total wavefunction is factorized into slow (nuclear) variables, which are treated classically, and fast ones (electrons) which remain quantum-mechanical. In the Ehrenfest method, nuclei move according to classical mechanics on a potential energy surface given by the expectation value of the electronic Hamiltonian. 

```math
i\hbar \dot{c}_i(t) = V_i(\mathbf{R}) c_i (t)
- i\hbar \sum_j \dot{\mathbf{R}} \cdot \mathbf{d}_{ij}(\mathbf{R})c_j(t)
```


## Example
Below the example of the Ehrenfest implementation is presented, using model from [Ananth2007](@cite).

First, the simulation parameters are created. Here we have a single atom with a mass of
`2000` and we are using Tully's third model, provided by [NonadiabaticModels.jl](@ref).

At the start, we assign `atoms` variable and initialise simulation using the atoms details and model employed by [NonadiabaticModels.jl](@ref).
```@example ehrenfest
using NQCDynamics

atoms = Atoms(1980)
sim = Simulation{Ehrenfest}(atoms, AnanthModelOne())
```
Next, we set dynamical distribution details:
```@example ehrenfest
e = 0.03
k = sqrt(e*2*atoms.masses[1])
r = Normal(-5, 1/sqrt(0.25))
v = k / atoms.masses[1]
distribution = DynamicalDistribution(v, r, size(sim1))* SingleState(1, Adiabatic())
```
To run an ensemble simulation we additionally choose number of trajectories `n_traj` and timespan `tspan` and we pass all the established settings to the `run_ensemble` function. In this example we output velocities by specifying `output=:velocity` and store the final values in the `final_velocities` array. Following that, we calculate final momenta.
```@example ehrenfest
n_traj = 500
tspan = (0.0, 3000.0)
solution1 = run_ensemble(sim1, tspan, distribution; 
    trajectories=n_traj, output=:velocity)
final_velocities = [r.velocity[end] for r in solution1]
momenta = reduce(vcat, final_velocities*atoms.masses[1])
```
Resulting momenta can be ploted by using `StatsPlots` package.
```@example ehrenfest
using StatsPlots
plot(density(momenta))
xlims!(-20,20)
```
