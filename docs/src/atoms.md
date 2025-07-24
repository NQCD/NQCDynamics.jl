```@setup logging
@info "Expanding src/atoms.md..."
start_time = time()
```
# [Handling Atoms](@id atoms)

!!! tip "Tip: NQCDynamics.jl handles atoms differently to `ase`"
	This package makes the choice to separate the atomic parameters from their positions and
	velocities for ease of use with the differential equations solvers.
	This contrasts somewhat with most other software packages where these would usually
	be joined together into a single object.

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
large vectors with many atom types.
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

When working with NQCDynamics, the most useful packages are [NQCBase](@ref)
for reading and writing structures and trajectories, and [PythonCall](https://github.com/JuliaPy/PythonCall.jl)
for working with [ASE](https://wiki.fysik.dtu.dk/ase/index.html) from within Julia. 

!!! tip "Tip: Installing Python dependencies"
	PythonCall.jl by default depends on CondaPkg.jl to install Python packages, but you can [opt out and use your own python installation](https://juliapy.github.io/PythonCall.jl/stable/pythoncall/#If-you-already-have-Python-and-required-Python-packages-installed).

## Using Python's `ase` package from Julia

Julia provides multiple options to run Python-based code from Julia scripts. [NQCDInterfASE](@ref) 
provides compatibility with [**PythonCall.jl**](https://github.com/JuliaPy/PythonCall.jl).
PythonCall forces you as a user to think more about when data is copied in memory between Python and Julia, supporting the creation of more efficient code. 

This example shows how `ase.build` can be used to build a structure from within Julia, then convert
from the ASE format into the required objects for atoms, positions and unit cell for an NQCD simulation:

```julia
using PythonCall
using NQCBase

ase_build = pyimport("ase.build")

# Make a silicon supercell using ASE
atoms_ase = ase_build.bulk("Si") * pytuple((4, 1, 1))
```

It is currently not possible to use an AtomsBase system directly with NQCDynamics, but can
be quickly converted to the correct format:

```julia
using NQCDynamics
using NQCDInterfASE # Import python interface functionality

atoms_nqcd, positions_nqcd, cell_nqcd = convert_from_ase_atoms(atoms_ase)

println(atoms_nqcd)
println(positions_nqcd)
println(cell_nqcd)
```


## Saving and loading with AtomsBase

After running a simulation it is often desirable to save the trajectory in a standard format for visualization.
For this, convert the NQCDynamics output into the AtomsBase format,
then use AtomsIO to write the file in your chosen format.

```julia
using AtomsIO
using NQCBase

v = rand(3, 8); #velocities
r = rand(3,8); #positions
cell = PeriodicCell([ 1 0 0;
         0 1 0;
         0 0 1;
            ]) # unit cell definition
#all definitions here are in atomic units, use Unitful.jl if you want to convert into Angstrom

system = System(atoms_nqcd, r, v, cell)

AtomsIO.save_system("Si.xyz", system) # Save a single image

trajectory = Trajectory(atoms_nqcd, [r, r, r, r], [v, v, v, v], cell)
AtomsIO.save_trajectory("Si.xyz", trajectory) # Save a trajectory
```

AtomsIO also provides `load_system` and `load_trajectory` which can be converted to the
NQCDynamics format as above to initialise simulations.
Refer to [AtomsIO](https://mfherbst.github.io/AtomsIO.jl/stable/) for more information.
```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
