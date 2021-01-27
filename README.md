# NonadiabaticMolecularDynamics.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maurergroup.github.io/NonadiabaticMolecularDynamics.jl/dev/)

Designing initial conditions for computer experiments

**TODO**
1. Diatomic state-to-state scattering (harmonic oscilator, rigid rotor, verticla position, kinetic energy) + thermalised surface (@Connor)
This would work for classical MD and MDEF (incl. RPMD extension)

2. Thermalised on-surface dynamics (molecule on surface positions/momenta sampled from Boltzmann distribution, long-time NVT simulation OR 
using West, Estreicher method: https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.96.115504 [already implemented in ASE_local by Reini])
This would work for classical MD and MDEF (incl. any ring-polymer extension) (@Reini)

3. Surface Hopping initial conditions(classical position/momenta sampling + select initial condition for electronic density matrix, and combine them)
This should work for the surface hopping case of 1 and 2

4. Monte Carlo sampling of initial conditions

5. State-selective analysis of trajectories for state-to-state scattering histogram binning (@Connor)


6. Make it possible to deliver initial conditions (momenta/positions) to different codes: ASE format, VENUS, Julia package?!

7. Excited state sampling: For small molecules we should do Wigner sampling, for large molecules classical MD and take every xth step. Every initial condition needs an excited state calulation. The oscillator strengths for every transition need to be calculated (energies and transition dipoles - this most often requires the reference method unless we have a model which can fit both values accurately) and according to the initial excitation wave length we set the initial state to start the dynamics. The initial state is the output which we then give to julia/MD via an input file(?). 
Worth to consider: If we want to do something with excited-states, maybe do the sampling in the beginning of a study and use these data points for fitting (at least when ML is applied).

![image](https://user-images.githubusercontent.com/38594562/100458375-91ea4880-30bb-11eb-9118-bb8eb97eaed3.png)
