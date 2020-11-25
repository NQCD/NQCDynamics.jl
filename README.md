# DICE.jl
Designing initial conditions for computer experiments

**TODO**
1. Diatomic state-to-state scattering (harmonic oscilator, rigid rotor, verticla position, kinetic energy) + thermalised surface (@Connor)
This would work for classical MD and MDEF (incl. RPMD extension)

2. Thermalised on-surface dynamics (molecule on surface positions/momenta sampled from Boltzmann distribution, long-time NVT simulation OR 
using West, Estreicher method: https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.96.115504 [already implemented in ASE_local by Reini])
This would work for classical MD and MDEF (incl. any ring-polymer extension)

3. Surface Hopping initial conditions(classical position/momenta sampling + select initial condition for electronic density matrix, and combine them)
This should work for the surface hopping case of 1 and 2

4. Monte Carlo sampling of initial conditions

5. State-selective analysis of trajectories for state-to-state scattering histogram binning (@Connor)


6. Make it possible to deliver initial conditions (momenta/positions) to different codes: ASE format, VENUS, Julia package?!
