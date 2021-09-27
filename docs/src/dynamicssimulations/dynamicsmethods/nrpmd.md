
# NRPMD

## Theory

Nonadiabatic ring polymer molecular dynamics is a method that uses the ring polymer
formalism to include quantum effects in the nuclear dynamics and mapping variables
for the electronic degrees of freedom.[^Richardson2013] [^Richardson2017] [^Chowdhury2019]
This results in a classical dynamics in an extended phasespace of the ring polymer
with each bead coupled to a set of classical mapping variables.
Originallly, this method was proposed as a simple combination of the 
Meyer-Miller-Stock-Thoss mapping formalism with RPMD but has since been
rigorously derived from nonadiabatic Matsubara dynamics.[^Chowdhury2021]

The classical Hamiltonian conserved by the NRPMD dynamics is given by

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

For mapping variable methods of this type, a symplectic algorithm [^Church2018]
exists which has the advantage of long time stability and can be easily combined
with the standard algorithms for ring polymer time-evolution.
For NRPMD we have implemented this algorithm using the Cayley modified ring polymer
propagator[^Korol2019] and obtain accurate and efficient dynamics.
For few beads, similar performance to the OrdinaryDiffEq.jl algorithms
is obtained, but as the number of beads increases this algorithm becomes more effective.

### Generating the initial distribution

Currently we provide only for nonequilibrium simulations where the nuclear part
of the distribution is separable from the electronic part.
Typically the nuclear distribution will be sampled using Langevin dynamics or Monte Carlo
sampling and the electronic variables are confined to a single electronic state.
This is appropriate for modelling photoexcitation dynamics but is not yet suitable
for equilibrium simulations.
This would require also sampling a thermal distribution for the mapping variables.

### Form of the Hamiltonian

Diabatic models defined in `NonadiabaticModels.jl` are of the appropriate form for
this method though they provide the potential as a single matrix, rather than separating
the state-dependent and independent parts.
It has been suggested that defining the Hamiltonian such that the trace of the diabatic
matrix is minimised leads to more accurate results, but here we have not done this for
simplicity.

## Example

Using NRPMD we can reproduce the figure 3a in the 2019 paper of Chowdhury
and Huo[^Chowdhury2019].

First we must generate a thermal ring polymer distribution in a harmonic potential.
A simple way to do this is to use Monte Carlo sampling for the positions and
obtain velocities from a Boltzmann distribution.

```@example nrpmd
using NonadiabaticMolecularDynamics

atom = Atoms(1)

sim = RingPolymerSimulation(atom, Harmonic(dofs=1), 4; temperature=1/16)

r0 = zeros(size(sim))
output = InitialConditions.ThermalMonteCarlo.run_advancedmh_sampling(sim, r0, 5e3, Dict(:X=>1.0))
velocities = InitialConditions.BoltzmannVelocityDistribution(1/16, masses(sim))

distribution = InitialConditions.DynamicalDistribution(velocities, output, size(sim);
                                     state=1, type=:diabatic)
```

We can check the distribution by plotting the phasespace diagram for each of the points
in our distribution:

```@example nrpmd
using CairoMakie

flat_position = vcat([p[:] for p in distribution.position]...)
flat_velocity = vcat([rand(distribution.velocity)[:] for p in 1:length(flat_position)]...)
scatter(flat_position, flat_velocity)
```

As usual, the simulation method is given as the type parameter `{NRPMD}` and
the simulation constructor is given the atoms, model, number of beads,
temperature and degrees of freedom.

```@example nrpmd
sim = RingPolymerSimulation{NRPMD}(atom, DoubleWell(γ=0.1), 4; temperature=1/16)
```

Next, we can use this distribution as a starting point for the dynamics simulations.
This will result in each trajectory starting from a random configuration in the
distribution.
Since the distribution also has `state=1` and `type=:diabatic`, the electronic
variables will be selected such that each trajectory is initialised in diabatic state 1.

The quantities output by the ensemble simulation are specified by the `Output` and
the `Reduction`.
The `Output` tells the simulation the quantities that we want to extract from each
trajectory, and the `Reduction` reduces the data.
By default the `Reduction` appends the output from each trajectory but here we would
like to average the results so we use the `MeanReduction`.

```@example nrpmd
output = Ensembles.OutputDiabaticPopulation(sim)
reduction = Ensembles.MeanReduction()
nothing # hide
```

Finally, we can combine the parameters and run the simulation.
The resulting plot shows the time dependent population difference and closely matches
the figure from the paper we were attempting to reproduce. Nice!

```@example nrpmd
ensemble = Ensembles.run_ensemble(sim, (0.0, 30.0), distribution; trajectories=1000,
                                  output=output, reduction=reduction, dt=0.1)

lines(0:0.1:30, [p[1]-p[2] for p in ensemble.u])
```

## References

[^Richardson2013]: Jeremy O. Richardson, Michael Thoss, [J. Chem. Phys. 139, 031102 (2013)](https://doi.org/10.1063/1.4816124)
[^Richardson2017]: Jeremy O. Richardson et al., [Chemical Physics 482 124–134 (2017)](https://doi.org/10.1016/j.chemphys.2016.09.036)
[^Chowdhury2019]: Sutirtha N. Chowdhury, Pengfei Huo, [J. Chem. Phys. 150, 244102 (2019)](https://doi.org/10.1063/1.5096276)
[^Chowdhury2021]: Sutirtha N. Chowdhury, Pengfei Huo, [J. Chem. Phys. 154, 124124 (2021)](https://doi.org/10.1063/5.0042136)
[^Church2018]: Matthew S. Church et al., [J. Chem. Phys. 148, 102326 (2018)](https://doi.org/10.1063/1.5005557)
[^Korol2019]: Roman Korol, Nawaf Bou-Rabee, Thomas F. Miller, [https://doi.org/10.1063/1.5120282](https://doi.org/10.1063/1.5120282)
