# [Fewest-switches surface hopping (FSSH)](@id fssh-dynamics)

Tully's FSSH [Tully1990](@cite) is one of the most popular methods for nonadiabatic
molecular dynamics and is classified as a mixed-quantum classical method, where the nuclei
are treated classically and the electrons are treated quantum mechanically.

The central concept that governs surface hopping methods is that the nuclei evolve on
a single adiabatic potential energy surface at any given moment.
At each timestep, a hopping probability is evaluated. If the hopping probability is larger than a uniform random number between 0 and 1,
the active state is switched and the adiabatic propagation continues on the new electronic state.
When this algorithm is applied to an ensemble of trajectories, the discrete adiabatic state 
populations approximate the quantum mechanical populations for each state.

The surface hopping classical Hamiltonian can be written as
```math
H(t) = \frac{1}{2} \mathbf{P}^T \mathbf{M}^{-1} \mathbf{P} + \sum_i \delta(s(t) - i) E_i(\mathbf{R})
```
where ``\mathbf{P}`` is the vector of momenta, ``\mathbf{R}`` the positions, and ``\mathbf{M}``
the diagonal mass matrix.
``s(t)`` can be viewed as a digital signal that takes on the value of the currently
occupied adiabatic state.
As such, this Hamiltonian describes classical dynamics that proceeds under the influence
of the potential ``E_i(\mathbf{R})`` when ``s(t) = i``.
The summation runs over all adiabatic states.

Of course, to integrate the associated equations of motion, ``s(t)`` must be obtained.
This quantity is obtained stochastically for each trajectory by making probabilistic hops
between surfaces.
The probabilities are obtained by integrating the electronic Schrödinger equation
alongside the dynamics as
```math
i\hbar \dot{c}_i(t) = E_i(\mathbf{R}) c_i (t)
- i\hbar \sum_j \dot{\mathbf{R}} \cdot \mathbf{d}_{ij}(\mathbf{R})c_j(t)
```
In this equation, ``c_i(t)`` are the complex coefficients for state ``i`` and
``\mathbf{d}_{ij}`` is the nonadiabatic coupling between adiabatic states ``i`` and ``j``.
The hopping probability is calculated as
```math
\gamma_{i \to j} = \sum_{\alpha} 2 \frac{P_\alpha}{M_\alpha}
\Re(\frac{\sigma_{ji}}{\sigma_{ii}}) d_{\alpha ij} dt.
```
At each timestep, a random number between 0 and 1 is generated which is compared to the
probabilities. If the probability is higher than the random number, then a hop is attempted.

Additionally in the fewest-switches scheme, the energy is conserved
for each trajectory by rescaling the momenta whenever a hop is performed.
As such, when a hop is attempted, it will only be successful when there is sufficient
kinetic energy for the energy to be conserved after the hop.
If there is insufficient kinetic energy, this is termed a frustrated hop, and
the dynamics proceeds without performing a hop.
When a hop is successful, the kinetic energy is adjusted and ``s(t)`` takes on the value
of the newly occupied state.
For a more detailed description of the algorithm and the momentum rescaling procedure, please
refer to [Subotnik2016](@cite). 
In this reference, the notion of reversing the momenta during frustrated hops is discussed.
In our implementation we leave the frustrated trajectories unchanged, though it is suggested
that the momentum reversal procedure may lead to better results in some cases.

## Algorithm

1. Integrate classical dynamics for one timestep
2. Integrate electronic dynamics for one timestep
3. Evaluate hopping probability
4. Perform hop if sufficient probability and kinetic energy
5. Rescale velocity if hop is performed
6. Return to step 1

!!! note

    With `DifferentialEquations.jl` we use a callback to perform the surface hopping
    procedure such that steps 1 and 2 are performed by the DifferentialEquations solvers and steps 3, 4, 5 are
    performed by the callback.

## Example

In this section we can investigate the results obtained for a single trajectory using FSSH.

First, the simulation parameters are created. Here, we have a single atom with a mass of
`2000` a.u. and we are using Tully's third model ([Tully1990](@cite)), provided by [NQCModels.jl](@ref).
```@example fssh
using Random; Random.seed!(10) # hide
using NQCDynamics

atoms = Atoms(2000)
sim = Simulation{FSSH}(atoms, TullyModelThree())
```

The [`DynamicsVariables`](@ref) constructor has some extra arguments for FSSH.
The first three match the classical case, but we also provide the initial state and
whether we want this state to be `Adiabatic()` or `Diabatic()`.
The type of state can be important when considering the ordering of the states.
The adiabatic states are always arranged from lowest to highest energy, whereas the diabatic
states will be ordered as defined in the model.
You can inspect the fields of `u` to ensure the initilisation has proceeded as you intend.
```@example fssh
u = DynamicsVariables(sim, [20/2000;;], [-10.;;], SingleState(1, Adiabatic()))
```

Finally, the trajectory can be run by passing all the parameters we have set up so far.
Here, we request both the discrete `:state` output which is equal to ``s(t)`` and 
`:population`, which gives us the population of each diabatic state along the trajectory.
```@example fssh
traj = run_trajectory(u, (0.0, 2000.0), sim, output=(:state, :population))
```

Now we can plot ``s(t)`` throughout the trajectory. The FSSH algorithm attempts to minimise
the total number of hops; in the limit of infinite hops the result would tend to the
mean-field (Ehrenfest) result, which is what FSSH attempts to avoid.
```@example fssh
using Plots
plot(traj, :state)
```

Similarly, we can plot the diabatic populations. Since FSSH is performed in the adiabatic
representation, even in the case of few hops, the diabatic populations can look dramatically
different depending on the complexity of the model Hamiltonian. 
```@example fssh
plot(traj, :population)
```

[Another example is available](@ref examples-tully-model-two) where we use FSSH and other
methods to reproduce some of the results from [Tully1990](@cite).
