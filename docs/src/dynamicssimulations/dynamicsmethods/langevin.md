# [Classical Langevin dynamics](@id langevin-dynamics)

Langevin dynamics can be used to sample the canonical ensemble for a classical system.
Langevin dynamics are based on classical equations of motion that are modified by
an additional drag force and a random force.
The Langevin equation of motion can be written as
```math
\mathbf{M}\ddot{\mathbf{R}} = - \nabla_R V(\mathbf{R}) + \mathbf{F}(t) - \gamma \dot{\mathbf{R}}
```
where ``\mathbf{M}`` are the masses, ``\ddot{\mathbf{R}}`` the time-derivative of the positions, 
i.e., the velocities, ``\nabla_R V(\mathbf{R})`` the gradient of the potential and ``\mathbf{F}(t)`` the random force
that is related to the friction coefficient ``\gamma`` by the second fluctuation-dissipation theorem.

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
As usual, ``\mathbf{P}`` is the vector of particle momenta and ``\mathbf{M}`` their diagonal mass matrix.

!!! note "Stochastic differential equations"

    There are two mathematical frameworks for handling stochastic differential equations,
    developed by Ruslan Stratonovich and Kiyosi Ito. To learn about the difference between the two
    in a physical context refer to [Risken1989](@cite).

As a stochastic differential equation, these two can be integrated immediately using
`StochasticDiffEq` provided by `DifferentialEquations`, which offers a variety of stochastic
solvers.
It is possible to exploit the dynamical structure of the differential equations
by splitting the integration steps into parts that can be solved exactly. In this context, 
it has been shown that the BAOAB method from [Leimkuhler2012](@cite) achieves good accuracy compared to other
similar algorithms and this algorithm is used here as the default.

## Example

Using Langevin dynamics we can sample the canonical ensemble for a simple harmonic
oscillator and investigate the energy expectation values.

Firstly we set up our system parameters. Here, we have two atoms in a harmonic
potential at a temperature of `1e-3`. We have arbitrarily chosen the dissipation constant
``\gamma = 1``, this can be tuned for optimal sampling in more complex systems. 
```@example langevin
using NQCDynamics
using Unitful

atoms = Atoms([:H, :C])
temperature = 1e-3
sim = Simulation{Langevin}(atoms, Harmonic(m=atoms.masses[1]); γ=1, temperature)
```

!!! note "Atomic units"

    As usual, all quantities are in atomic units by default.

Here we can generate a simple starting configuration with zeros for every degree of freedom.
```@example langevin
u = DynamicsVariables(sim, zeros(size(sim)), zeros(size(sim)))
```

Running the dynamics proceeds by providing all the parameters along with
any extra keywords. This time we have requested both the positions and velocities to be
outputted and have selected a timestep `dt`.
Since the default algorithm is a fixed timestep algorithm an error will be thrown if a
timestep is not provided.
```@example langevin
traj = run_dynamics(sim, (0.0, 2000.0), u; output=(OutputPosition, OutputVelocity), dt=0.1)
```

Here, we plot the positions of our two atoms throughout the simulation.
```@example langevin
using Plots
plot(traj, :OutputPosition, label=["Hydrogen" "Carbon"], legend=true)
```

We next plot the velocities. Notice how the carbon atom with its heavier mass has a smaller
magnitude throughout.
```@example langevin
plot(traj, :OutputVelocity, label=["Hydrogen" "Carbon"], legend=true)
```

Using the configurations from the Langevin simulation we can obtain expectation values along
the trajectories.
This can be done manually, but we provide the [`Estimators`](@ref) module to make this
as simple as possible.

!!! note Estimators

    [Here](@ref `Estimators`) you can find the available quantities that [`Estimators`](@ref) provides.
    To add new quantities, you must implement a new function inside `src/Estimators.jl`.

Let's find the expectation for the potential energy during our simulation.
This is the potential energy of the final configuration in the simulation:
```@repl langevin
Estimators.potential_energy(sim, traj[:OutputPosition][end])
```
We could evaluate this for every configuration and average it manually.
Fortunately however, we have the [`@estimate`](@ref Estimators.@estimate) macro that
will do this for us:
```@repl langevin
Estimators.@estimate potential_energy(sim, traj[:OutputPosition])
```

!!! tip

    We can verify this result by comparing to the
    equipartition theorem which states that each quadratic degree of freedom should contribute
    ``\frac{1}{2}kT`` to the total energy.
    As this is a harmonic system, this gives us the exact classical potential energy
    expectation as equal to the temperature, since we have two degrees of freedom and
    we are in atomic units.

Similarly, we can evaluate the kinetic energy expectation with:
```@repl langevin
Estimators.@estimate kinetic_energy(sim, traj[:OutputVelocity])
```
Again, this takes a similar value since the total energy is evenly split between the kinetic
and potential for a classical harmonic system.
