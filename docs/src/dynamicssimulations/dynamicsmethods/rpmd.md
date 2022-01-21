# [Ring polymer molecular dynamics (RPMD)](@id rpmd-dynamics)

Ring polymer molecular dynamics is a quantum dynamics methods that attempts
to approximate Kubo-transformed real-time correlation functions ([Craig2004](@cite)).

The idea is to exploit the classical isomorphism that maps a quantum particle onto
the extended phasespace of a classical ring polymer.
By this we mean that it can be shown that the quantum partition function
for a system can be manipulated such that it resembles the classical partition function
of a system containing many replicas of the original particle joined to together with
harmonic springs in a ring.
In the limit of infinite beads or replicas in the ring polymer, the isomorphism becomes
exact and it is possible to evaluate quantum expectation values by evaluating
ensemble averages for the classical ring polymer system.
This is referred to as the field of imaginary-time path integrals and the techniques used
are Path Integral Monte Carlo (PIMC) and Path Integral Molecular Dynamics (PIMD)
depending on whether molecular dynamics or Monte Carlo methods are used to explore the
phasespace.

RPMD was proposed as a heuristic extension of imaginary-time path integrals to evaluate
real-time dynamical quantities.
To perform RPMD, it is necessary to solve Hamilton's equations for the ring polymer
Hamiltonian:
```math
H = \sum_\alpha^N \frac{1}{2} \mathbf{P}_\alpha^T \mathbf{M} \mathbf{P}_\alpha
+ \frac{1}{2} \omega_N^2
(\mathbf{R}_\alpha - \mathbf{R}_{\alpha+1})^T
\mathbf{M}
(\mathbf{R}_\alpha - \mathbf{R}_{\alpha+1})
+ V(\mathbf{R}_\alpha)
```
where the ring polymer spring constant ``\omega_N = 1 / \hbar\beta_N`` and
``\beta_N = \beta / N``.

When the initial distribution is taken as the thermal ring polymer distribution and
this Hamiltonian is used to generate configurations at later times,
the correlation functions obtained can be used to approximate real-time quantum correlation
functions.

## Example

Let us perform some simple adiabatic ring polymer dynamics to get a feel
for what the ring polymer dynamics looks like.

```@example rpmd
using NonadiabaticMolecularDynamics
using Unitful

atoms = Atoms([:H])
sim = RingPolymerSimulation(atoms, Free(2), 50; temperature=100u"K")
```

We have set up a 2D system by giving the [`Free`](@ref) model 2 degrees of freedom and
have specified that the ring polymer should have 50 beads.

We can initialise the simulation with zero velocity and a random distribution for the
ring polymer bead positions. For a real RPMD simulation you will use the thermal ring
polymer distribution obtained from a PIMC or Langevin simulation but here for simplicity
we use a normally distributed configuration.
```@example rpmd
u = DynamicsVariables(sim, zeros(size(sim)), randn(size(sim)))
nothing # hide
```

Now we can run the simulation:
```@example rpmd
dt = 2.5
traj = run_trajectory(u, (0.0, 500.0), sim; output=(:position), dt=dt)
nothing # hide
```

We can visualise this ring polymer trajectory with a 2D scatter plot that shows how
the ring polymer evolves in time. Here we have joined the adjacent beads together with
lines, with the end and start beads joined with a different color.
Hopefully you can see that this interaction appears equivalent to all the others.

```@example rpmd
using CairoMakie

rs = traj.position

index = Node(1)
xs = @lift(rs[$index][1,1,:])
ys = @lift(rs[$index][2,1,:])
close_loop_x = @lift([rs[$index][1,1,end], rs[$index][1,1,begin]])
close_loop_y = @lift([rs[$index][2,1,end], rs[$index][2,1,begin]])
fig = scatter(xs, ys, axis = (title = @lift("t = $(round(Int, dt*($index-1)))"),))
lines!(xs, ys)
lines!(close_loop_x, close_loop_y)
xlims!(-3, 3)
ylims!(-3, 3)

timestamps = 1:length(traj.position)
record(fig, "../../assets/figures/rpmd.gif", timestamps;
        framerate = 30) do i
    index[] = i
end
nothing
```

![rpmd fig](../../assets/figures/rpmd.gif)

!!! note

    We have used Makie's animation features to produce this animation. If you want
    information how this stuff works, take a look at the
    [Makie documentation](https://makie.juliaplots.org/stable/documentation/animation/).

Since this package is focused on nonadiabatic dynamics, you won't see much adiabatic RPMD
elsewhere in the documentation but it's useful to understand how the original adiabatic
version works before moving onto the nonadiabatic extensions.
