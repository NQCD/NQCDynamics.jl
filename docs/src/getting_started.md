# Getting Started

Here you can find a brief overview of how to run a simple dynamics simulation.

## Simple classical dynamics

Let's create the atoms in our simulation.
We can see that `Atoms` contains the symbols, atomic numbers and masses for our two atoms.
```@example classical
using NonadiabaticMolecularDynamics # hide
atoms = Atoms([:H, :C])
```

Next we should define the potential energy surface our forces will be derived from.
This is a simple harmonic potential that imposes a harmonic force on each our atoms.
```@example classical
model = Harmonic(ω=50.0, r₀=0.0)
```

We can now combine these two objects inside our `Simulation`.
Additionally, here we have specified that each atom has only a single degree of freedom.
```@example classical
sim = Simulation{Classical}(atoms, model; DoFs=1)
nothing # hide
```

The final step before solving our trajectory is to specify the initial velocities
and positions for each atom.
```@example classical
z = DynamicsVariables(sim, randn(1, 2), randn(1, 2))
```

Finally, let us run the trajectory.
The second argument here is the time span for the dynamics.
By default, classical dynamics uses the `VelocityVerlet` algorithm and we must provide a
timestep.
```@example classical
solution = Dynamics.run_trajectory(z, (0.0, 50.0), sim, dt=0.1)
nothing # hide
```

`Dynamics.run_trajectory` returns a `Table` from `TypedTables.jl`.
This contains the time and the values of the output quantities at each time.
By default, it outputs the value of the `DynamicalVariables` into the field `u`.

Convenience functions are provided that allow easy plotting of the result:
```@example classical
using Plots
plot(solution, :u)
```

## Moving forward
So far we have walked through a simple example that illustrates the solution
of a single trajectory for a model system.
Later sections of the documentation detail the specifics associated with different dynamics
methods but each will closely follow the above format.
