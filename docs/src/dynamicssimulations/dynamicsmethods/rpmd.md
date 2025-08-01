```@setup logging
@info "Expanding src/dynamicssimulations/dynamicsmethods/rpmd.md..."
start_time = time()
```
# [Ring polymer molecular dynamics (RPMD)](@id rpmd-dynamics)

Ring polymer molecular dynamics is a quantum dynamics methods that attempts
to approximate Kubo-transformed real-time correlation functions ([Craig2004](@cite)).

The idea is to exploit the classical isomorphism that maps a quantum particle onto
the extended phasespace of a classical ring polymer.
It can be shown that the quantum partition function
for a system can be manipulated such that it resembles the classical partition function
of a system containing many replicas of the original particle joined to together with
harmonic springs in a ring.
In the limit of infinite beads or replicas in the ring polymer, the isomorphism becomes
exact and it is possible to evaluate quantum expectation values by evaluating
ensemble averages for the classical ring polymer system.
This is referred to as the field of imaginary-time path integrals and the techniques used
are Path Integral Monte Carlo (PIMC) and Path Integral Molecular Dynamics (PIMD)
depending on whether molecular dynamics or Monte Carlo methods are used to explore the
phasespace ([tuckerman2010](@cite)).

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
We set up a 2D system for one hydrogen atom by giving the [`Free`](@ref) model 2 degrees of freedom and
specify that the ring polymer should have 50 beads.

```@example rpmd
using NQCDynamics
using Unitful

atoms = Atoms([:H])
sim = RingPolymerSimulation(atoms, Free(2), 50; temperature=100u"K")
```

!!! note "Atomic units"

    Recall that the quantities are always in atomic units unless [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/)
    has been used to specify alternative units. The temperature here has been specified using Kelvin.

We initialise the simulation with zero velocity and a random distribution for the
ring polymer bead positions. For a real RPMD simulation you will use the thermal ring
polymer distribution obtained from a PIMC or Langevin simulation but here for simplicity
we use a normally distributed configuration.
```@example rpmd
u = DynamicsVariables(sim, zeros(size(sim)), randn(size(sim)))
nothing # hide
```

!!! tip

    To learn how to work with the thermal ring polymer phase space, refer to the [Storing and sampling distributions](@ref nqcdistributions) section.

Now we can run the simulation, for which we use the time interval 0.0 to 500.0 and a time 
step of `dt = 2.5`:
```@example rpmd
dt = 2.5
traj = run_dynamics(sim, (0.0, 500.0), u; output=OutputPosition, dt=dt)
nothing # hide
```

We can visualise this ring polymer trajectory with a 2D scatter plot that shows how
the ring polymer evolves in time. Here, we have joined the adjacent beads together with
lines, with the end and start beads joined with a different color.
This animation shows the cyclic nature of the ring polymer, and how every bead is connected
to its two neighbours.

```@example rpmd
using Plots

rs = traj[:OutputPosition]

timestamps = 1:length(traj[:OutputPosition])
@gif for i in timestamps
    xs = rs[i][1,1,:]
    ys = rs[i][2,1,:]
    close_loop_x = [rs[i][1,1,end], rs[i][1,1,begin]]
    close_loop_y = [rs[i][2,1,end], rs[i][2,1,begin]]

    plot(
        xlims=(-3, 3),
        ylims=(-3, 3),
        legend=false
    )

    plot!(xs, ys, color=:black)
    scatter!(xs, ys)
    plot!(close_loop_x, close_loop_y)
end
```

Since this package is focused on nonadiabatic dynamics, you won't see much adiabatic RPMD
elsewhere in the documentation, but it's useful to understand how the original adiabatic
version works before moving on to the nonadiabatic extensions.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
