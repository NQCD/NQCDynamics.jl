# [General Introduction](@id mdef-dynamics)
 
A set of fundamental and technologically relevant chemical processes (surface scattering, dissociative chemisorption, surface diffusion, recombinative desorption, etc) are often catalyzed at the metal surface of several late transition metals (Au, Ag, Cu, Pt, Pd, Rh, etc). These metallic surfaces, unlike to other kind of surfaces, are characterized by highly dense electronic states landscape which produce a continue conduction and valence bands virtually without any band gap. A proper theoretical description of these chemical processes are often challenging due to the Born-Oppenheimer (BO) approximation is not longer valid and non-adiabatic effects have to be considered to describe the energy exchange that can take place between adsorbate and substrate degrees of freedom (DOF).

A fully quantum dynamic approach of this complex scenario is currently unfeasible and the gas-surfaces reaction dynamics are often approached by quasi-classical methods where the nuclear motion are treated as classical particles but a proper electronic structure description of the metal surface is included through first principle electronic structure calculations in different ways. 

## Molecular dynamics with electronic friction (MDEF)

Molecular dynamics with electronic friction (MDEF) is one of main workhorse used to deal with the non-adiabaticity in gas-surfaces chemical reactions and it has been widely employed to decribe and simulate the nuclear dynamics in several molecular systems. MDEF is a theoretical model based on ground-state generalized Langevin equation (GLE) of motion which allow to introduce non-adiabatic effects by means of both friction and stochastic forces. This approach was originally introduced by Head-Gordon and Tully and the non-adiabatic effects can be included through different electronic friction models (see section below, LDFA and TDPT). Within this theoretical framework, the coupling of molecular degree of freedom to electron-hole pair (EHPs) excitations in the metal are described by means of a frictional force  which condense the metal substrate electronic structure into electronic friction object (a tensor or a coefficient). In the context of GLE, the nuclear coordinates of each adsorbate atom evolves as follow, 

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \}) }{\partial \mathbf{r_{i}}}  -{f_{r,i}^{fric}}  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```
The first term on the right hand side of the equation (1) corresponds to conservative force associated with potential energy surface (PES) as in the adiabatic case. The second term is the friction force and it come from multiplication between the electronic friction object (``f_{e,i}^{fric}``) and the velocity. This means that the final friction force contribution depends dramatically of both quantities. Finally, the last term is temperature and friction-dependent stochastic random force which assure the detailed balance. For some particular cases, the random force can be neglected setting the electronic temperature at 0 K (see scattering event example below).

### Reactive scattering events: Initial conditions & compuatational settings

The current version of our Julia implementation allow us to simulate the state-to-state vibrational de-excitation probability on reactive scattering events at metal surfaces. To run this kind of simulations, a set of initial positions and velocities with a specific set of ro-vibrational quantum states ``\nu`` and ``j`` have to be initially generated (see initial condition section/ sampling methods/ semiclassical EBK quantisation). With a specific ro-vibrational quantum state is possible to compute different properties after molecular collision and energy transfer with the metal surface like the vibrational de-excitation probabilities. A full description of our reactive system also includes other basic variables definition like type of atoms, unit cell and electronic temperature which can be properly defined before to run any simulation.

```@docs
system = load("distribution.jld2")

Tel = 0u"K"
cell = system["cell"]
distribution = system["dist"]
atoms = system["atoms"]

```
Here, "distribution.jld2" is a binary file which was previously generated (see initial conditions section) and contains relevant information about simulated system like the initial distribution (positions and velocities), substrate unit cell and types of adsorbate atoms. Each one of these variables can be accessed providing the "key" string to get the assocaited "value" due to dictionary structure associted with jld2 files. Also, the electronic temperature (``T_el``) can be selected for our MDEF simulations. For the reactive scattering event the electronic temperature was set equal to 0 negleting the last term (random force) in the equation 1.  

The set of velocities and positions stored on "distribution.jld2" ("system ["dist"]") can be read in two different ways, OrderedSelection or RandomSelection. The first one reads the velocities and positons in a ordered way starting from the first member of the list and Randomselction randomly selects the positions and velocities. In our examples, we have used the orderedSelections but in other situation the RandomSelection can be useful.

```@docs
selection = Ensembles.OrderedSelection(distribution)

```
As other molecular dynamics simulations, the simulated total time (tspan) and time steps (dt) can be suitably selected depending on the specific conditons and simulated system, 

```@docs
tspan = (0.0, 420.0u"fs")
dt=0.1u"fs"

```
Finally, depending on the simulated system is possible to include a "termination function" which stops the simulation when a set of conditions are satified. This feature can be especially useful for scattering events to reduce the computational cost and prevent any unphysical situation.


```@docs
function termination_condition(u, t, integrator)::Bool
   add_conditions
end

terminate = Dynamics.TerminatingCallback(termination_condition)

```

## Local density friction approximation (LDFA)

Local density friction approximation (LDFA) is a theoretical model which describes the electronic friction ``(f_{e,i}^{fric})`` term in the above equation based on the electron density of the metal substrate. This approximation assumes a single friction coefficient (``\eta_{e,i}``) for each adsorbate atom assuming a anisotropic contribution. In the LDFA theoretical framework the above equation read as

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \})}{\partial \mathbf{r_{i}}}   -\eta_{e,i}(\mathbf{r_{i}})  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```

In our current LDFA implementation, a set of pre-calculated electronic friction coefficients (``\eta_{e,i}``) computed at different Wigner-Seitz radius (``r_s``) are used to fit and get an analytical expression to connect any ``r_s`` values with an single electronic friction coefficient by means of   
cubic Spline functions. The Wigner-Sietz radius is connected with the metal substrate electron density by the following equation, 

```math
   r_s(\rho) = (\frac{3}{4\pi \rho (\mathbf{r_{i}})})^{1/3}
```

In this way, the electron density associated with the current substrate atom position is used to compute the respective friction coefficient through fitting function for each point of the trajectory.

### Computational settings

To run any MDEF simulation, we need first to provide a model to compute both the electronic friction elements and the conservative forces from a potential surface energy (PES). In the case of LDFA simulation is also necesary to provide a "cube" file which is used to get the elecronic friction coefficients by using the fittng function.


```@docs
 model = H2AgModel()
 model=LDFAModel(model,"density.cube",atoms,friction_atoms=[1,2],cell)

```
Here, the models variable contains the potential, force and friction data. For this example, "H2AgModel" contains the data associated with PES and friction machine learning models pre-calculated. The second line only modify the way the friction object is computed in the context of LDFA model. Also, it is necesarry to selected the atoms which we want to compute the frictions coefficients.

Now, we have all the necesary parameters to run our MDEF simulations, 

```@docs
sim = Simulation{MDEF}(atoms, model, cell=cell,temperature=Tel)
ensemble = Ensembles.run_ensemble(sim,tspan,selection;dt=0.1u"fs",trajectories=75000,output=Ensembles.OutputFinal(),callback=terminate,ensemble_algorithm=EnsembleDistributed())
```
 

## Time-dependent Perturbation theory (TDPT)

A more general formulation of the electronic friction object was also developed under the umbrella of electronic friction tensor(EFT) or orbital-dependent electronic friction (ODF) approaches. Both formulations are essentially equivalent and they incorporate the anisotropy nature of the electronic friction object by a multidimentional tensor (``\Lambda_{ij}``) instead of a single coefficient as usually computed at LDFA level.  The electronic friction elements can be computed by first-principle calculations in the context of first-order time-dependent perturbation theory (TDPT) at the density functional theory (DFT) level. In the context of this approach the GLE read as, 

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \}) }{\partial \mathbf{r_{i}}}   -\sum_{j} \Lambda_{ij}  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```

Each electronic friction tensor (EFT) elements corresponds to relaxation rate due to electron-nuclear coupling along the Cartesian coordinate ``i`` due to motion in the ``j`` direction. The electronic friction tensor elements can be computed by using the Fermi's golden rule.

Here, ``\vert \psi_{k\nu} \rangle`` and ``\epsilon_{k\nu}`` are the Kohn-Sham (KS) ground state eigenstates and eigenenergies, respectively. The derivatives quatities are computed by finite difference numerically and normalized Gaussian distribution of finite width (``\sigma``) centered at Fermi level is used to facilitate convergence instead to the ``\delta`` function. A ``\delta`` value of 0.6 is often selected to due is able to produce converged results in the majority of the systems analyzed.

``\Lambda_{ij}`` is object with (``3N\times3N``)-dimension where N is often the total number of adsorbate atoms considered explicitly on the study system.

### Computational settings

In this case a simular setting than LDFA is used but without the cube file. 

```@docs
 model = H2AgModel()

sim = Simulation{MDEF}(atoms, model, cell=cell,temperature=Tel)
ensemble = Ensembles.run_ensemble(sim,tspan,selection;dt=0.1u"fs",trajectories=75000,output=Ensembles.OutputFinal(),callback=terminate,ensemble_algorithm=EnsembleDistributed())
```


```@autodocs
Modules=[Dynamics]
Pages=["Dynamics/mdef.jl"]
```
