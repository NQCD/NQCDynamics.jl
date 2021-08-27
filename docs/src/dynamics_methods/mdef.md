# General Introduction
A set of fundamental and technologically relevant chemical processes are often activated at the metal surface of several late transition metals (Au, Ag, Cu, Pt, Pd, Rh, etc.) Metal surfaces, unlike to other surfaces, are characterized by highly dense electronic state landscape which produce a continue conduction and valence bands virtually without any band gap. The gas-surface reactions dynamics at metal surfaces often involves a set of different chemical events (scattering, dissociative chemisorption, surface diffusion, recombinative desorption, etc) that characterize the considered system and the intrinsic catalytic and electrochemical properties associated with the analyzed metal substrate. A proper theoretical description of this set of processes is often challenging due to the Born-Oppenheimer (BO) approximation is not longer valid and non-adiabatic effects have to be considered to describe the potential energy exchange that can take place between adsorbate and substrate degrees of freedom (DOF).

A fully quantum dynamic approach of this complex scenario is currently unfeasible and this kind of process (gas-surfaces reaction dynamis) is often approached by quasi-classical methods where the nuclear motion are treated as classical particles but a proper electronic structure description of the metal is included from first principle electronic structure calculations. 

# MDEF
Molecular dynamics with electronic frinction (MDEF) is one of main workhorse used to deal with the non-adiabaticity in gas-surfaces chemical reaction and it has been widely employed to decribe the nuclear dynamics in several surface process. MDEF is a theoretical model based on ground-state generalized Langevin equation (GLE) of motion which allow to introduce non-adiabatic effect by means of friction and stochastic forces. This approach was originally introduced by Head-Gordon and Tully  and improves the initial description based on Newton equation to describe gas-surface reaction at metal surface by including non-adiabatic effects through different electronic friction models. In this theoretical framework, the coupling of molecular degree of freedom to electron-hole pair (EHPs) excitation in the metal are described by means of a frictional force  which condense the metal substrate electronic structure into electron friction component.  In the context of GLE, the temporal evolution of the nuclear dynamics for each adsorbate atom is not only governed by the potential energy surface (V) but also for two extra terms as is shown below


# LDFA

# TDPT

```@autodocs
Modules=[Dynamics]
Pages=["Dynamics/mdef.jl"]
```
