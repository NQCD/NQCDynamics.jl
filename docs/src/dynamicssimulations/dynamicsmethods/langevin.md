# [Classical Langevin dynamics](@id langevin-dynamics)

Langevin dynamics can be used to sample the canonical ensemble for a classical system.
It involves a classical dynamics that is modified with an additional drag force
and a random force.
The equation of motion can be written as
```math
\mathbf{M}\ddot{\mathbf{R}} = - \nabla_R V(\mathbf{R}) + \mathbf{F}(t) - \gamma \dot{\mathbf{R}}
```
where the friction constant ``\gamma`` is related to the vector of white noise processes
``\mathbf{F}(t)`` by the second fluctuation-dissipation theorem.

Equally the above equation can be written in the form of Ito stochastic differential
equations [Leimkuhler2012](@cite)
```math
\begin{aligned}
d\mathbf{R} &= \mathbf{M}^{-1} \mathbf{P} dt\\
d\mathbf{P} &= [-\nabla V(\mathbf{R}) - \gamma \mathbf{P}] dt
+ \sigma \mathbf{M}^{1/2} d\mathbf{W}
\end{aligned}
```
where ``\sigma = \sqrt{2\gamma/\beta}`` and ``\mathbf{W}`` is a vector of ``N`` independent
Wiener processes.

As a stochastic differential equation, this can be integrated immediately using
`StochasticDiffEq` provided by `DifferentialEquations` using a variety of stochastic
solvers.
However, it is possible to exploit the dynamical structure of the differential equations
by splitting the integration steps into parts that can be solved exactly.
It has been shown previously that the BAOAB method achieves good accuracy compared to other
similar algorithms and this algorithm is used here as the default. ([Leimkuhler2012](@cite))

## Example

Using Langevin dynamics we can sample the canonical ensemble for a simple harmonic
oscillator and investigate the energy expectation values.

Firstly we should set up our system parameters. here we have two atoms in a harmonic
potential at temperature of `1e-3`. We have arbitrarily chosen the dissipation constant
``\gamma = 1``, this can be tuned for optimal sampling in more complex systems. 
```@example langevin
using NQCDynamics
using Unitful

atoms = Atoms([:H, :C])
sim = Simulation{Langevin}(atoms, Harmonic(m=atoms.masses[1]); γ=1, temperature=1e-3)
```

Here we can generate a simple starting configuration with zeros for every degree of freedom.
```@example langevin
u = DynamicsVariables(sim, zeros(size(sim)), zeros(size(sim)))
```

Running the dynamics proceeds in the usual way by providing all the parameters along with
any extra keywords. This time we have requested both the positions and velocities to be
output and have selected a timestep `dt`.
Since the default algorithm is a fixed timestep algorithm an error will be thrown if a
timestep is not provided.
```@example langevin
traj = run_trajectory(u, (0.0, 2000.0), sim; output=(:position, :velocity), dt=0.1)
```

Here we can see the positions of our two atoms throughout the simulation.
```@example langevin
using Plots
plot(traj, :position)
```

These are the velocities, notice how the carbon atom with its heavier mass has a smaller
magnitude throughout.
```@example langevin
plot(traj, :velocity)
```

With these configurations we can calculate expectation values by computing averages along
the trajectories.
This can be done manually, but we provide the [`Estimators`](@ref) module to make this
as simple as possible.

Let's find the expectation for the potential energy during our simulation.
This is the potential energy of the final configuration in the simulation:
```@repl langevin
Estimators.potential_energy(sim, traj.position[end])
```
We could evaluate this for every configuration and average it manually.
Fortunately however, we have the [`@estimate`](@ref Estimators.@estimate) macro that
will do this for us:
```@repl langevin
Estimators.@estimate potential_energy(sim, traj.position)
```

!!! tip

    We can verify this result by comparing to the
    equipartition theorem which states that each quadratic degree of freedom should contribute
    ``\frac{1}{2}kT`` to the total energy.
    Since this is a harmonic system, this gives us the exact classical potential energy
    expectation as equal to the temperature, since we have two degrees of freedom and
    we are in atomic units.

Similarly, we can evaluate the kinetic energy expectation with:
```@repl langevin
Estimators.@estimate kinetic_energy(sim, traj.velocity)
```
Again, this takes a similar value since the total energy is evenly split between the kinetic
and potential for a classical harmonic system.
