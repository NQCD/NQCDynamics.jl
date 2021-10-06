# [General Introduction](@id mdef-dynamics)
 
A set of fundamental and technologically relevant chemical processes (surface scattering, dissociative chemisorption, surface diffusion, recombinative desorption, etc) are often catalyzed at the metal surface of several late transition metals (Au, Ag, Cu, Pt, Pd, Rh, etc). These metallic surfaces, unlike to other kind of surfaces, are characterized by highly dense electronic state landscape which produce a continue conduction and valence bands virtually without any band gap. A proper theoretical description of these chemical processes are often challenging due to the Born-Oppenheimer (BO) approximation is not longer valid and non-adiabatic effects have to be considered to describe the energy exchange that can take place between adsorbate and substrate degrees of freedom (DOF).

A fully quantum dynamic approach of this complex scenario is currently unfeasible and the gas-surfaces reaction dynamics are often approached by quasi-classical methods where the nuclear motion are treated as classical particles but a proper electronic structure description of the metal is included through first principle electronic structure calculations in different ways. 

# Molecular dynamics with electronic friction (MDEF)

MDEF is one of main workhorse used to deal with the non-adiabaticity in gas-surfaces chemical reactions and it has been widely employed to decribe and simulate the nuclear dynamics in several systems. MDEF is a theoretical model based on ground-state generalized Langevin equation (GLE) of motion which allow to introduce non-adiabatic effects by means of both friction and stochastic forces. This approach was originally introduced by Head-Gordon and Tully and the non-adiabatic effects can be included through different electronic friction models (see section below). Within this theoretical framework, the coupling of molecular degree of freedom to electron-hole pair (EHPs) excitation in the metal are described by means of a frictional force  which condense the metal substrate electronic structure into electronic friction object. In the context of GLE, the nuclear coordinates for each adsorbate atom evolves as follow, 

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \}) }{\partial \mathbf{r_{i}}}  -{f_{r,i}^{fric}}  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```
The first term on the right hand side of the equation (1) corresponds to conservative force associated with potential energy surface (PES) like the adiabatic case. The second term is friction force and it come from multiplication between the electronic friction object (``f_{e,i}^{fric}``) and the velocity. This means that the final force contribution depends dramatically of both quantities. Finally, the last term is temperature and friction-dependent stochastic random force which assure the detailed balance. For some particular cases, the random force can be neglected setting the electronic temperature at 0 K.

### Reactive scattering events: Initial conditions

State-to-state vibrational de-excitation probability on reactive scattering events at metal surface can be properly simulated within the ecosystem of our implementation. To run this kind of simulations, a set of initial positions and velocities with a specific set of ro-vibrational quantum states ``\nu`` and ``j`` should be initially generated (see Initial Condition section/sampling methods/Semiclassical EBK quantisation).  Likewise, other basic variables like atoms, cell and electronic temperature should be defined before to run any simulation.

```@docs
system = load("distribution.jld2")

Tel = 0u"K"
cell = system["cell"]
distribution = system["dist"]
atoms = system["atoms"]

```
Here, "distribution.jld2" binary file was previously generated (see initial conditions) and contains relevant information about simulated system like the initial distribution (positions and velocities), substrate unit cell and types of adsorbate atoms. Each one of these variables can be accessed providing the "key" to get the assocaited values due to dictionary structure associted with jld2 files. Also, the electronic temperature (Tel) can be selected for our MDEF simulations. 

The set of velocities and positions stored on "distribution.jld2" (system["dis"]) can be read in two different ways, OrderedSelection or RandomSelection. The first one reads the velocities and positons in a ordered way starting from the first members of the list and Randomselction selects randomly the positions and velocities. In our examples, we have used the orderedSelections.

```@docs
selection = Ensembles.OrderedSelection(distribution)

```
As other MD simulations, the total time and time steps can be suitably selected depending on the specific conditons used. 

```@docs
tspan = (0.0, 420.0u"fs")
dt=0.1u"fs"

```
Finally, our implementation allow us to include a termination function which stops the simulation when a set of conditions are satified. This feature can be especially useful for scattering events to reduce the computational cost. 


```@docs
function termination_condition(u, t, integrator)::Bool
   add_conditions
end

terminate = Dynamics.TerminatingCallback(termination_condition)

```

# Local Density Friction Approximation (LDFA)

LDFA is a theoretical model which describes the electronic friction ``(f_{e,i}^{fric})`` term in the above equation based on the electron density of the metal substrate. This approximation assumes a single friction coefficient (``\eta_{e,i}``) for each adsorbate atom assuming a anisotropic contribution. In the LDFA theoretical framework the above equation read as

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \})}{\partial \mathbf{r_{i}}}   -\eta_{e,i}(\mathbf{r_{i}})  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```

In our current LDFA implementation, a set of pre-calculated electronic friction coefficients (``\eta_{e,i}``) computed at different Wigner-Seitz radius (``r_s``) are used to fit and get an analytical expression to connect any ``r_s`` values with an single electronic friction coefficient by means of   
cubic Spline functions. The Wigner-Sietz radius is connected with the metal substrate electron density by the following equation, 

```math
   r_s(\rho) = (\frac{3}{4\pi \rho (\mathbf{r_{i}})})^{1/3}
```

In this way, the electron density associated with the current substrate atom position is used to compute the respective friction coefficient through fitting function for each point of the trajectory.

### Settings

To run a LDFA simulation is necesary to provide a "cube" file for the substrate which is used to get the elecronic friction coefficients. For any MDEF simulation a model and type of simulations need to be declared.

```@docs
 model = H2AgModel()
 model=LDFAModel(model,"density.cube",atoms,friction_atoms=[1,2],cell)

```
Here, models contains potential, force and friction data. For this example, "H2AgModel" contains the data associated with PES and friction machine learning models pre-calculated. The second line only modify the way the friction object is computed in the context of LDFA model. Also, it is necesarry to selected the atoms which we want to compute the frictions coefficients.

Now, we have all the necesary parameters to run our LDFA simulations, 

```@docs
sim = Simulation{MDEF}(atoms, model, cell=cell,temperature=Tel)
ensemble = Ensembles.run_ensemble(sim,tspan,selection;dt=0.1u"fs",trajectories=75000,output=Ensembles.OutputFinal(),callback=terminate,ensemble_algorithm=EnsembleDistributed())
```


### 

# Time-dependent Perturbation theory (TDPT)

A more general formulation of the electronic friction object was also developed under the umbrella of \textit{electronic friction tensor}(EFT) or \textit{orbital-dependent electronic friction} (ODF) approaches. Both formulations are essentially equivalent and they incorporate the anisotropy nature of the electronic friction object by a multidimentional tensor (``\Lambda_{ij}``) instead of a single coefficient as usually computed at LDFA level.  The electronic friction elements can be computed by first-principle calculations in the context of \textit{first-order time-dependent perturbation theory} (TDPT) at the density functional theory (DFT) level. In the context of this approach the GLE read as, 

```math
   m_{i} \frac{d^{2} \mathbf{r_{i}} }{dt^{2}} = -\frac{\partial V (\{ \mathbf{r_{j}} \}) }{\partial \mathbf{r_{i}}}   -\sum_{j} \Lambda_{ij}  \frac{d \mathbf{r_{i}} }{dt} + \mathbf{{R_{i}}}(t)
```

Each electronic friction tensor (EFT) elements corresponds to relaxation rate due to electron-nuclear coupling along the Cartesian coordinate ``i`` due to motion in the ``j`` direction. The electronic friction tensor elements can be computed by using the Fermi's golden rule.\cite{2016f_maurer,2016_maurer}

Here, ``\vert \psi_{k\nu} \rangle`` and ``\epsilon_{k\nu}`` are the Kohn-Sham (KS) ground state eigenstates and eigenenergies, respectively.\cite{2016f_maurer,2016_maurer} The derivatives quatities are computed by finite difference numerically and normalized Gaussian distribution of finite width (``\sigma``) centered at Fermi level is used to facilitate convergence instead to the ``\delta`` function. A ``\delta`` value of 0.6 is often selected to due is able to produce converged results in the majority of the systems analyzed.

``\Lambda_{ij}`` is object with (``3N\times3N``)-dimension where N is often the total number of adsorbate atoms considered explicitly on the study system.

### Settings

To run a LDFA simulation is necesary to provide a "cube" file for the substrate which is used to get the elecronic friction coefficients. For any MDEF simulation a model and type of simulations need to be declared.

```@docs
 model = H2AgModel()

```
Here, models contains potential, force and friction data. For this example, "H2AgModel" contains the data associated with PES and friction machine learning models pre-calculated. The second line only modify the way the friction object is computed in the context of LDFA model. Also, it is necesarry to selected the atoms which we want to compute the frictions coefficients.

Now, we have all the necesary parameters to run our LDFA simulations, 

```@docs
sim = Simulation{MDEF}(atoms, model, cell=cell,temperature=Tel)
ensemble = Ensembles.run_ensemble(sim,tspan,selection;dt=0.1u"fs",trajectories=75000,output=Ensembles.OutputFinal(),callback=terminate,ensemble_algorithm=EnsembleDistributed())
```




```@autodocs
Modules=[Dynamics]
Pages=["Dynamics/mdef.jl"]
```
