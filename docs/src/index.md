# NonadiabaticMolecularDynamics.jl

This is a package for performing nonadiabatic molecular dynamics simulations, along with
generating the initial conditions and analysing the resulting trajectories.

Dynamics methods:

- Classical molecular dynamics
- Classical Langevin dynamics
- Fewest-switches surface hopping
- Molecular dynamics with electronic friction
- Ring polymer molecular dynamics

Initial conditions sampling:

- Metropolis-Hastings Markov chain Monte Carlo thermal sampling

An advantage of this package is that analytic model Hamiltonians and atomic systems are
treated equivalently within the dynamics, making the transition from models to
atomic systems as frictionless as possible.
