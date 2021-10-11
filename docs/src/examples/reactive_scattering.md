# Reactive scattering from a metal surface

The current version of our Julia implementation allow us to simulate the state-to-state
vibrational de-excitation probability on reactive scattering events at metal surfaces.
To run this kind of simulations, a set of initial positions and velocities with a specific
set of ro-vibrational quantum states ``\nu`` and ``j`` have to be initially generated
(see initial condition section/ sampling methods/ semiclassical EBK quantisation).
With a specific ro-vibrational quantum state is possible to compute different properties
after molecular collision and energy transfer with the metal surface like the vibrational
de-excitation probabilities.
A full description of our reactive system also includes other basic variables definition
like type of atoms, unit cell and electronic temperature which can be properly defined
before to run any simulation.

```@example mdef
system = load("distribution.jld2")

Tel = 0u"K"
cell = system["cell"]
distribution = system["dist"]
atoms = system["atoms"]

```
Here, "distribution.jld2" is a binary file which was previously generated
(see initial conditions section) and contains relevant information about simulated system
like the initial distribution (positions and velocities), substrate unit cell and types of
adsorbate atoms.
Each one of these variables can be accessed providing the "key" string to get the associated
 "value" due to dictionary structure associted with jld2 files.
 Also, the electronic temperature (``T_el``) can be selected for our MDEF simulations.
 For the reactive scattering event the electronic temperature was set equal to 0 neglecting
 the last term (random force) in the equation 1.  

The set of velocities and positions stored on "distribution.jld2" ("system ["dist"]") can
be read in two different ways, OrderedSelection or RandomSelection.
The first one reads the velocities and positons in a ordered way starting from the first
member of the list and Randomselction randomly selects the positions and velocities.
In our examples, we have used the orderedSelections but in other situation the
RandomSelection can be useful.

```@example mdef
selection = Ensembles.OrderedSelection(distribution)

```
As other molecular dynamics simulations, the simulated total time (tspan) and time steps
(dt) can be suitably selected depending on the specific conditons and simulated system, 

```@example mdef
tspan = (0.0, 420.0u"fs")
dt=0.1u"fs"

```
Finally, depending on the simulated system is possible to include a "termination function" which stops the simulation when a set of conditions are satified. This feature can be especially useful for scattering events to reduce the computational cost and prevent any unphysical situation.


```@example mdef
function termination_condition(u, t, integrator)::Bool
   add_conditions
end

terminate = Dynamics.TerminatingCallback(termination_condition)

```

### Computational settings

To run any MDEF simulation, we need first to provide a model to compute both the electronic friction elements and the conservative forces from a potential surface energy (PES). In the case of LDFA simulation is also necesary to provide a "cube" file which is used to get the elecronic friction coefficients by using the fittng function.


```@example mdef
 model = H2AgModel()
 model=LDFAModel(model,"density.cube",atoms,friction_atoms=[1,2],cell)

```
Here, the models variable contains the potential, force and friction data. For this example, "H2AgModel" contains the data associated with PES and friction machine learning models pre-calculated. The second line only modify the way the friction object is computed in the context of LDFA model. Also, it is necesarry to selected the atoms which we want to compute the frictions coefficients.

Now, we have all the necesary parameters to run our MDEF simulations, 

```@example mdef
sim = Simulation{MDEF}(atoms, model, cell=cell,temperature=Tel)
ensemble = Ensembles.run_ensemble(sim,tspan,selection;dt=0.1u"fs",trajectories=75000,output=Ensembles.OutputFinal(),callback=terminate,ensemble_algorithm=EnsembleDistributed())
```
 

### Computational settings

In this case a similar setting than LDFA is used but without the cube file. 

```@example mdef
 model = H2AgModel()

sim = Simulation{MDEF}(atoms, model, cell=cell,temperature=Tel)
ensemble = Ensembles.run_ensemble(sim,tspan,selection;dt=0.1u"fs",trajectories=75000,output=Ensembles.OutputFinal(),callback=terminate,ensemble_algorithm=EnsembleDistributed())
```

