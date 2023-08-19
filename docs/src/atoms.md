# [Handling Atoms](@id atoms)

This package makes the choice to separate the atomic parameters from their positions and
velocities for ease of use with the differential equations solvers.
This contrasts somewhat with most other software packages where these would be usually
by joined together into a single object.

The atomic parameters here are contained within the
[`Atoms`](@ref) type introduced earlier
in the [Getting started](@ref) section.
As mentioned previously, there exist some basic constructors which use either elemental
symbols or numbers to initialise the parameters:
```@repl atoms
using NQCDynamics
Atoms([:H, :H, :H])
Atoms([1, 2, 3])
```

If there are many atoms, you can use Julia's array manipulation utilities to create
large vectors with many atoms types.
For example, if adding an adsorbate to a metal surface, it could be initialised as:
```@repl atoms
au = fill(:Au, 40)
no = [:N, :O]
auno = [au; no]
Atoms(auno)
```

# [Manipulating atomic structures with AtomsBase.jl](@id atoms-base)

[AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) provides a convenient format for
representing atomic geometries, facilitating interoperability between a collection of
different packages.

When working with NQCDynamics, the most useful packages are [AtomsIO](https://github.com/mfherbst/AtomsIO.jl)
for reading and writing structures and trajectories, and [ASEconvert](https://github.com/mfherbst/ASEconvert.jl)
for working with [ASE](https://wiki.fysik.dtu.dk/ase/index.html) from within Julia.

## Using ASE with ASEconvert.jl

This example shows how ASEconvert can be used to build a structure, then convert
from the ASE format into an AtomsBase compatible system:

```@example atomsbase
using ASEconvert

# Make a silicon supercell using ASE
atoms_ase = ase.build.bulk("Si") * pytuple((4, 1, 1))

# Convert to an AtomsBase-compatible structure
atoms_ab = pyconvert(AbstractSystem, atoms_ase)
```

It is currently not possible to use an AtomsBase system directly with NQCDynamics, but can
be quickly converted to the correct format:

```@repl atomsbase
using NQCDynamics
atoms_nqcd = Atoms(atoms_ab)
r = Position(atoms_ab)
v = Velocity(atoms_ab)
c = Cell(atoms_ab)
```

## Saving and loading with AtomsIO.jl

After running a simulation it often desirable to save the trajectory in a standard format for visualization.
For this, convert the NQCDynamics output into the AtomsBase format,
then use AtomsIO to write the file in your chosen format.

```@example atomsbase
using AtomsIO

system = System(atoms_nqcd, r, v, c)

AtomsIO.save_system("Si.xyz", system) # Save a single image

trajectory = Trajectory(atoms_nqcd, [r, r, r, r], [v, v, v, v], c)
AtomsIO.save_trajectory("Si.xyz", trajectory) # Save a trajectory
```

AtomsIO also provides `load_system` and `load_trajectory` which can be converted to the
NQCDynamics format as above to initialize simulations.
Refer to [AtomsIO](https://mfherbst.github.io/AtomsIO.jl/stable/) for more information.
