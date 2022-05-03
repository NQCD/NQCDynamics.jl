# [Nonadiabatic ring polymer molecular dynamics (NRPMD)](@id nrpmd-dynamics)

## Theory

Nonadiabatic ring polymer molecular dynamics (NRPMD) is a method that uses the ring polymer
formalism to include quantum effects in the nuclear dynamics and mapping variables
for the electronic degrees of freedom.
([Richardson2013](@cite), [Richardson2017](@cite), [Chowdhury2019](@cite))
This results in a classical dynamics in an extended phasespace of the ring polymer
with each bead coupled to a set of classical mapping variables.
Originally, this method was proposed as a simple combination of the 
Meyer-Miller-Stock-Thoss mapping formalism with RPMD but has since been
rigorously derived from nonadiabatic Matsubara dynamics ([Chowdhury2021](@cite)).

The classical Hamiltonian conserved by NRPMD is given by

```math
H_N = \sum_{\alpha=1}^N \left[
\frac{P_\alpha^2}{2M} + V_0(R_\alpha)
+ \frac{M}{2\beta_N^2\hbar^2} (R_\alpha - R_{\alpha-1})^2
+ \frac{1}{\hbar} \sum_{nm}V_{nm}(R_\alpha)
([q_\alpha]_n[q_\alpha]_m + [p_\alpha]_n[p_\alpha]_m - \delta_{nm}\hbar)
\right]
```

which contains ``N`` replicas with positions ``R_\alpha`` and momenta ``P_\alpha``
joined by harmonic springs.
``M`` is the mass, ``\beta_N = \beta / N`` is the inverse temperature scaled by the number of beads.
``V_0(R_\alpha)`` is the state independent potential.
``V_{nm}`` are the matrix elements of the diabatic potential. The sum runs over all pairs of states.
Each replica has a set of mapping variables ``[q_\alpha]_n`` and ``[p_\alpha]_n``
that interact only within the set associated with a single replica.
The consequence of this is that the electronic dynamics is not contaminated by interbead 
coupling.

The equations of motion obtained from this Hamiltonian are

```math
\begin{aligned}
\dot{R}_\alpha &= \frac{P_\alpha}{M}
\\
\dot{P}_\alpha &=
- \frac{M}{\beta_N^2 \hbar^2}(2 R_\alpha - R_{\alpha+1} - R_{\alpha-1})
- \nabla_{R_\alpha} V_0(R_\alpha)
- \frac{1}{2\hbar}\sum_{nm} \nabla_{R_\alpha} V_{nm}(R_\alpha)
([q_\alpha]_n[q_\alpha]_m + [p_\alpha]_n[p_\alpha]_m - \delta_{nm}\hbar)
\\
[\dot{q}_\alpha]_n &=
\frac{1}{\hbar} \sum_m V_{nm}(R_\alpha)[p_\alpha]_m
\\
[\dot{p}_\alpha]_n &=
-\frac{1}{\hbar} \sum_m V_{nm}(R_\alpha)[q_\alpha]_m
\end{aligned}
```

## Implementation details

### Solving the differential equations

For mapping variable methods of this type, a symplectic algorithm ([Church2018](@cite))
exists. This algorithm has the advantage of long time stability and can be easily combined
with the standard algorithms for ring polymer time-evolution.
For NRPMD we have implemented this algorithm using the Cayley modified ring polymer
propagator ([Korol2019](@cite)) and obtain accurate and efficient dynamics.
For few beads, similar performance to the OrdinaryDiffEq.jl algorithms
is obtained, but as the number of beads increases this algorithm becomes more effective.

### Generating the initial distribution

Currently, we provide this functionality only for nonequilibrium simulations where the nuclear part
of the distribution is separable from the electronic part.
Typically, the nuclear distribution will be sampled using Langevin dynamics or Monte Carlo
sampling and the electronic variables are confined to a single electronic state.
This is appropriate for modelling photoexcitation dynamics but is not yet suitable
for equilibrium simulations.
Equilibrium dynamics would require also sampling a thermal distribution for the mapping variables.

### Form of the Hamiltonian

The diabatic models defined in `NQCModels.jl` are of the appropriate form for
this method though they provide the potential as a single matrix, rather than separating
the state-dependent and independent parts.
It has been suggested that defining the Hamiltonian such that the lowest eigenvalue
of the diabatic matrix is zero everywhere leads to improved convergence in the sampling
([Richardson2013](@cite)).
However, here we have not done this for simplicity when defining the models.

## Example

Using NRPMD we can reproduce the Fig. 3a in the 2019 paper of Chowdhury
and Huo ([Chowdhury2019](@cite)).

First we generate a thermal ring polymer distribution in a harmonic potential.
A simple way to do this is to use Monte Carlo sampling for the positions and
obtain velocities from a Boltzmann distribution.

```@example nrpmd
using NQCDynamics

atom = Atoms(1)

sim = RingPolymerSimulation(atom, Harmonic(dofs=1), 4; temperature=1/16)

r0 = zeros(size(sim))
steps = 5e3 # Number of Monte Carlo steps
step_size = Dict(:X=>1.0) # Monte Carlo step size for species :X
output = InitialConditions.ThermalMonteCarlo.run_advancedmh_sampling(sim, r0, steps, step_size)
velocities = VelocityBoltzmann(1/16, masses(sim), (1,1))

distribution = DynamicalDistribution(velocities, output, size(sim)) * PureState(1)
```

!!! note "`size(sim)`"

    `size(sim)` returns the system size as `(degrees of freedom, number of atoms, number of beads)`.

!!! tip "Monte Carlo sampling"

    Further information on Monte Carlo sampling can be found [here](@ref mhmc-sampling).

We can check the distribution by plotting the phasespace diagram for each of the points
in our distribution:

```@example nrpmd
using CairoMakie

nuclear = distribution.nuclear
flat_position = reduce(vcat, (nuclear.position[i][:] for i in 1:length(nuclear)))
flat_velocity = reduce(vcat, (rand(nuclear.velocity)[:] for _ in 1:length(nuclear)))
scatter(flat_position, flat_velocity)
```

!!! note "`reduce(vcat, ...)`"

    Here we have used `reduce` in combination with `vcat` to vertically concatenate all of the information
    into a single array for plotting.

The simulation method is given as the type parameter `{NRPMD}` and
the simulation constructor is given the atoms, model, number of beads,
temperature and degrees of freedom.

```@example nrpmd
sim = RingPolymerSimulation{NRPMD}(atom, DoubleWell(γ=0.1), 4; temperature=1/16)
```

Next, we can use this distribution as a starting point for the dynamics simulations.
This will result in each trajectory starting from a random configuration in the
distribution.
For NRPMD, the electronic variables are sampled from a Gaussian, independent of the
initial electronic state.
The electronic state is introduced in the correlation function expression when correlating
the initial and final populations.

The quantities output by the ensemble simulation are specified by the `output` and the
`reduction`.
The `output` follows the `DifferentialEquations` format where we provide a function
that determines the output of each trajectory.
The `reduction` can be one of `:mean`, `:append`, or `:sum`, which will determine
how the data from each trajectory is combined.

!!! note "Ensemble simulations"

    Further details on ensemble simulations are available [here](@ref ensembles).

```@example nrpmd
output = TimeCorrelationFunctions.PopulationCorrelationFunction(sim, Diabatic())
nothing # hide
```

Finally, we can combine the parameters and run the simulation.
The resulting plot shows the time dependent population difference and closely matches
the figure from the paper we were attempting to reproduce. Nice!

```@example nrpmd
ensemble = run_ensemble(sim, (0.0, 30.0), distribution;
    trajectories=1000, output, reduction=:mean, dt=0.1,
    u_init=[zeros(2,2) for i=1:length(0:0.1:30.0)])

plt = lines(0:0.1:30, [p[1,1]-p[2,1] for p in ensemble])
plt.axis.xlabel = "Time"
plt.axis.ylabel = "Population difference"
plt
```
