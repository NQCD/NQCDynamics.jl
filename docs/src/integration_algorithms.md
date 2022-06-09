# Integration algorithms

At the core of NQCDynamics.jl is the DifferentialEquations.jl package that performs all of the dynamics simulations.
Within the sub-packages OrdinaryDiffEq.jl and StochasticDiffEq.jl, a variety of integration algorithms have been implemented
that are available to use without needing to implement custom algorithms for specific applications. 
However, in some cases it can be desirable to implement algorithms that can take advantage of the special structure of the dynamical system at hand.
A key example in the field of molecular dynamics is the famous velocity Verlet algorithm that is extremely popular due to its requirement for only
a single force evaluation during each time step and symplectic energy conservation properties.
In fact, velocity Verlet, along with a variety of other symplectic solvers are also implemented within OrdinaryDiffEq.jl.

For some problems encountered within semiclassical adiabatic and nonadiabatic dynamics,
there are a few different algorithms that can be used to obtain improved performance.
Ideally, these would be implemented using the DifferentialEquations.jl interface to allow for others to easily use these algorithms for their own problems.
However, it can be challenging to implement algorithms with an appropriate level of generality.
NQCDynamics.jl contains implementations for a few algorithms using the DifferentialEquations.jl interface
but the implementations are not completely generic and are coupled to the rest of the package.
In particular, the algorithms rely on specific formats for the dynamical arrays and use some functions that are not provided within the `DEProblem`.
In future it would be great to try to decouple the algorithms and package them separately within the DifferentialEquations.jl ecosystem so that others can use them more easily.

This page describes a few applications where special algorithms are available to enhance performance.
In each section it is noted which algorithms are available within NQCDynamics.jl.

## Ring polymer propagation

Path integral molecular dynamics and ring polymer molecular dynamics involve solving Hamilton's equations for a classical ring polymer Hamiltonian.
The ring polymer Hamiltonian describes many replicas of the system joined together by harmonic springs.
The number of replicas or beads required must be increased until convergence is achieved.
After adding many beads the ring polymer dynamics becomes hard to integrate, as the ring polymer internal modes become the highest
frequency modes in the system, limiting the largest acceptable time step.

To circumvent this difficulty, the ring polymer equations of motion can be partitioned to separate
the free ring polymer dynamics from the influence of the external potential.
Since the free ring polymer dynamics is entirely harmonic, it is possible to solve this part analytically,
allowing for time steps that are not limited by the internal ring polymer frequencies.
Recently, the symplectic Cayley modified algorithm has been demonstrated to exhibit strong stability and outperform the original algorithm.

!!! note

    The Cayley modified algorithm was originally introduced by [Korol2019](@cite).
    This paper provides a detailed description of ring polymer dynamics, how the integration algorithm works, 
    and benchmarks the performance of the algorithms. 

The Hamiltonian ring polymer integration algorithm has also been extended for thermostatted dynamics:
such as for thermal sampling in path integral molecular dynamics, or in thermostatted ring polymer molecular dynamics.
As an extension to the Cayley modified algorithm for the Hamiltonian dynamics, the work of [Korol2020](@cite) suggests
the BCOCB algorithm as the most effective for ring polymer dynamics with Langevin thermostatting.
The BCOCB nomenclature refers to the sub-steps within each time step.
B is the external potential, C is the Cayley modified free ring polymer step, and O is the thermostat.
In this nomenclature, the integration algorithm in the absence of the thermostat can be referred to as the BCB algorithm.

!!! tip

    NQCDynamics.jl implements both the BCOCB and BCB algorithms for Langevin and Hamiltonian dynamics, respectively.
    They are the default algorithms when performing adiabatic ring polymer dynamics.

## Mixed quantum-classical propagation

Mixed quantum-classical methods such as mean-field Ehrenfest dynamics or surface hopping dynamics involve simultaneous
propagation of nuclear and electronic sub-systems.
The two sub-systems evolve on different timescales and it can be advantageous to use different time steps or algorithms for each part.
Commonly the nuclear part is solved using the velocity Verlet algorithm and the electronic part is handled using a Runge-Kutta method or an exponential integrator.

!!! note 

    SHARC and Newton-X, two popular surface hopping codes, use the split-algorithm approach.
    In SHARC, the nuclear degrees of freedom are propagated using the velocity Verlet algorithm,
    whilst the wavefunction is propagated using an exponential integrator with a smaller time step. 
    Since an exponential integrator is exact when the propagation operator is constant,
    reducing the time step would have no benefit if the nuclei remained fixed during the time step.
    SHARC instead linearly interpolates the propagation operator during the electronic steps so that the nuclei propagation operator changes smoothly during the nuclear time step.
    This procedure is explained in the SHARC manual.
    Newton-X allows a few choices for the wavefunction integration algorithm and the nuclear quantities are interpolated similarly to SHARC.

NQCDynamics.jl uses the standard library of OrdinaryDiffEq.jl solvers to run the dynamics for mixed quantum-classical methods.
For model Hamiltonians where the evaluation of the electronic quantities is fast it is not necessary to use an augmented Verlet algorithm,
instead it is sufficient to use any of the adaptive solvers from OrdinaryDiffEq.jl.
However, in future it would be useful to implement some of these partioned algorithms that are able to achieve enhanced performance for large, expensive systems.

For ring polymer mixed quantum-classical methods, it is possible to combine the algorithms used for ring polymer propagation with the partitioning idea from the mixed quantum-classical solvers.
Even for model systems, the performance is significantly improved when the ring polymer modes are solved separately such that larger time steps can be used.
NQCDynamics.jl implements an augmented form of the BCB algorithm that uses the BCB algorithm for the ring polymer degrees of freedom and uses the Tsit5 algorithm from OrdinaryDiffEq.jl for the electronic part.
Currently the time steps for both sub-systems are fixed to be the same, but in future this constraint should be removed.
In fact, it would even be possible to use an adaptive solver for the electronic part that can automatically adjust the time step as necessary.

## Semiclassical mapping variable propagation

Mapping variable methods describe the coupled nuclear-electronic problem using a classical Hamiltonian,
where additional variables have been introduced to represent the electronic populations.
As is the case with the ring polymer methods,
it is possible to construct symplectic algorithms where each timestep is partitioned into exactly soluble sub-steps.
Of particular note is the MInt algorithm described by [Church2018](@cite).
The algorithm is symplectic, symmetric and time-reversible and can also be combined with the ring polymer algorithms for ring polymer dynamics with mapping variables.
NQCDynamics.jl implements the MInt algorithm and a variant of the MInt algorithm for ring polymer systems that combines the BCB algorithm for the nuclei with the MInt algorithm for the mapping variables.

## Electronic friction propagation

Molecular dynamics with electronic friction is described by Langevin equations of motion,
equivalent to classical Hamiltonian dynamics with an additional drag force and stochastic force.
In the case of Langevin dynamics with a constant friction coefficient, there has been much interest in developing algorithms of low order that can be used for large molecular dynamics simulations.
The work of [Leimkuhler2012](@cite) has demonstrated that, although a few choices exist, the BAOAB algorithm performs most favourably.

!!! tip

    The BAOAB algorithm uses the same nomenclature as the ring polymer algorithms introduced in [Ring polymer propagation](@ref).
    Each letter represents one component of a single time step. B is the external force step that updates the velocities, A is the position update, and O is the thermostatting step.
    For ring polymer dynamics the A step encompasses the free ring polymer step.

For electronic friction dynamics, the friction is described by a tensor, not a single number as for traditional Langevin dynamics.
However, the BAOAB algorithm is still applicable, but the O step requires the matrix exponential of the friction tensor.
Since the tensor is positive semi-definite, it is possible to perform the exponentiation by first diagonalising the tensor.
Using this approach, NQCDynamics.jl implements the BAOAB algorithm for tensorial friction and the BCOCB when using a ring polymer system. 

!!! note

    The `DynamicalSDEProblem` in StochasticDiffEq.jl represents systems that contain positions and velocities and have a stochastic component.
    The `DynamicalSDEProblem` was originally implemented for performing Langevin thermostatted dynamics simulations using the BAOAB algorithm.
    At the time of writing, BAOAB is the only algorithm implemented in StochasticDiffEq.jl for these problems.
    In future it would be useful to implement further algorithms and allow for more general noise profiles.
