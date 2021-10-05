# General Introduction

A set of fundamental and technologically relevant chemical processes (scattering, dissociative chemisorption, surface diffusion, recombinative desorption, etc) are often catalyzed at the metal surface of several late transition metals (Au, Ag, Cu, Pt, Pd, Rh, etc). These metallic surfaces, unlike to other kind of surfaces, are characterized by highly dense electronic state landscape which produce a continue conduction and valence bands virtually without any band gap. A proper theoretical description of this set of chemical processes are often challenging due to the Born-Oppenheimer (BO) approximation is not longer valid and non-adiabatic effects have to be considered to describe the energy exchange that can take place between adsorbate and substrate degrees of freedom (DOF).

A fully quantum dynamic approach of this complex scenario is currently unfeasible and this kind of processes (gas-surfaces reaction dynamics) are often approached by quasi-classical methods where the nuclear motion are treated as classical particles but a proper electronic structure description of the metal is included through first principle electronic structure calculations in different ways. 

# Molecular dynamics with electronic friction (MDEF)

MDEF is one of main workhorse used to deal with the non-adiabaticity in gas-surfaces chemical reactions and it has been widely employed to decribe the nuclear dynamics in several surface process. MDEF is a theoretical model based on ground-state generalized Langevin equation (GLE) of motion which allow to introduce non-adiabatic effect by means of both friction and stochastic forces. This approach was originally introduced by Head-Gordon and Tully and the non-adiabatic effects can be included through different electronic friction models (see section below). In this theoretical framework, the coupling of molecular degree of freedom to electron-hole pair (EHPs) excitation in the metal are described by means of a frictional force  which condense the metal substrate electronic structure into electronic friction object.  In the context of GLE, the temporal evolution of the nuclear coordinates for each adsorbate atom is not only governed by the potential energy surface (V) but also for two extra terms as is shown below, 

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \}) }{\partial \mathbf{r_{i}}}  -{f_{r,i}^{fric}}  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
\label{gle}
```
The first term on the right hand side of the equation (1) corresponds to conservative force associated with potential energy surface (PES) like the adiabatic case. The second term is friction force and it come from multiplication between the electronic friction object ($f_{e,i}^{fric}$) and the velocity. This means that the final force contribution depends dramatically of both quantities. Finally, the last term is temperature and friction-dependent stochastic random force which assure the detailed balance. For some particular cases, the random force can be neglected setting the electronic temperature at 0 K.

example,

# Local Density Friction Approximation 

LDFA is a theoretical model which describes the electronic friction (($f_{e,i}^{fric}$)) term in the above equation based on the electron density of the metal substrate. This approximation assumes a single friction coefficient ($\eta_{e,i}$) for each adsorbate atom assuming a anisotropic contribution. In the LDFA theoretical framework the above equation read as

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \})}{\partial \mathbf{r_{i}}}   -\eta_{e,i}(\mathbf{r_{i}})  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```

In our current LDFA implementation, a set of pre-calculated electronic friction coefficients ($\eta_{e,i}$) computed at different Wigner-Seitz radius ($r_s$) are used to fit and get an analytical expression to connect any $r_s$ values with an single electronic friction coefficient by means of   
cubic Spline functions. The Wigner-Sietz radius is connected with the metal substrate electron density by the following equation, 

```math
   r_s(\rho) = (\frac{3}{4\pi \rho (\mathbf{r_{i}})})^{1/3}
```

In this way, the electron density associated with the current substrate atom position is used to compute the respective friction coefficient through fitting function for each point of the trajectory.

example, 

# Time-dependent Perturbation theory

A more general formulation of the electronic friction object was also developed under the umbrella of \textit{electronic friction tensor}(EFT) or \textit{orbital-dependent electronic friction} (ODF) approaches. Both formulations are essentially equivalent and they incorporate the anisotropy nature of the electronic friction object by a multidimentional tensor ($\Lambda_{ij}$) instead of a single coefficient as usually computed at LDFA level.  The electronic friction elements can be computed by first-principle calculations in the context of \textit{first-order time-dependent perturbation theory} (TDPT) at the density functional theory (DFT) level. In the context of this approach the GLE read as, 

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \}) }{\partial \mathbf{r_{i}}}   -\sum_{j} \Lambda_{ij}  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```

Each \textbf{electronic friction tensor}(EFT) elements corresponds to relaxation rate due to electron-nuclear coupling along the Cartesian coordinate $i$ due to motion in the $j$ direction. The electronic friction tensor elements can be computed by using the Fermi's golden rule.\cite{2016f_maurer,2016_maurer}

Here, $\vert \psi_{k\nu} \rangle$ and $\epsilon_{k\nu}$ are the Kohn-Sham (KS) ground state eigenstates and eigenenergies, respectively.\cite{2016f_maurer,2016_maurer} The derivatives quatities are computed by finite difference numerically and normalized Gaussian distribution of finite width ($\sigma$) centered at Fermi level is used to facilitate convergence instead to the $\delta$ function. A $\delta$ value of 0.6 is often selected to due is able to produce converged results in the majority of the systems analyzed.

$\Lambda_{ij}$ is object with ($3N\times3N$)-dimension where N is often the total number of adsorbate atoms considered explicitly on the study system.

example,




```@autodocs
Modules=[Dynamics]
Pages=["Dynamics/mdef.jl"]
```
