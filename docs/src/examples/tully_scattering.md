# [Scattering probabilities for TullyModelTwo](@id examples-tully-model-two)

In this section we can roughly reproduce the results of figure 5 from [Tully1990](@cite).
This figure presents the scattering outcomes when a particle interacts with Tully's model 2
with an increasing magnitude of incident kinetic energy.

First, let's set up our system parameters:
```@example fssh
using NonadiabaticMolecularDynamics

atoms = Atoms(2000)
sim = Simulation{FSSH}(atoms, TullyModelTwo())
```

Each data point in the figure is obtained from an ensemble average of ...
