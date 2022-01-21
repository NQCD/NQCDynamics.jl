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

## [Reading and writing atomic structures](@id reading-and-writing)

When using a complex system however, it is likely more effective to read structures directly
from a file.
We provide two ways to do this, either using the 
[ExtXYZ.jl](https://github.com/libAtoms/ExtXYZ.jl) package which works for xyz files.
Or instead there is a conversion to and from the `ase.Atoms` type which can be used
when [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) is loaded.

First we can use `ase` to build a system and write it to an `.xyz` file.
```@example atoms
using PyCall

build = pyimport("ase.build")

slab = build.fcc100("Al", size=(2, 2, 3))
build.add_adsorbate(slab, "Au", 1.7, "hollow")
slab.center(axis=2, vacuum=4.0)

slab.write("slab.xyz")
```

Now we can read it in with the [`read_extxyz`](@ref NQCBase.read_extxyz)
function.
```@repl atoms
atoms, positions, cell = read_extxyz("slab.xyz")
atoms
positions
cell
```

Similarly, we can write the file with
[`write_extxyz`](@ref NQCBase.write_extxyz):
```@repl atoms
write_extxyz("out.xyz", atoms, positions, cell)
```
Both of these functions also work with trajectories such that the positions will be a vector
of configurations, whilst the atoms and cell will remain unchanged.

If not using `.xyz` files, we can directly use the IO capability of ase to read or the write
the files.
This can be done by using the conversions between our data types and the `ase.Atoms` object.
```@repl atoms
atoms = Atoms([:H, :H, :C])
ase_atoms = NQCBase.convert_to_ase_atoms(atoms, rand(3, 3))
NQCBase.convert_from_ase_atoms(ase_atoms)
```
These conversions work both ways such that you can read any file format using
ase then convert the `ase.Atoms` object to our types afterwards.
Then at the end when you are finished, you can convert them back and write your output
with ase.
