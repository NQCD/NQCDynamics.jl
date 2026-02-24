```@setup logging
@info "Expanding src/NQCModels/fullsizemodels.md..."
start_time = time()
```

# Full dimensional model library

## ASE interface

The easiest way to obtain potentials and forces from established codes is to
use the interfaces implemented in [ASE](https://wiki.fysik.dtu.dk/ase/).

We provide the [`ClassicalASEModel`](@ref) which wraps an ASE atoms object and its
associated calculator to implement the required [`potential`](@ref) and
[`derivative`](@ref) functions.
Several examples for connecting common machine-learning interatomic potentials to NQCModels.jl through the ASE interface are shown in the [MLIP examples](@ref ml-pes-models) section.

!!! note

    The interface works by calling the relevant Python functions using
    [PythonCall](https://github.com/JuliaPy/PythonCall.jl).
    To use PythonCall, you must make sure that your `python` version contains all the
    relevant packages, such as `ase`.
    [PythonCall](https://github.com/JuliaPy/PythonCall.jl) can be configured to use a particular
    pre-installed Python or install its own.
    Refer to the [PythonCall README](https://github.com/JuliaPy/PythonCall.jl) for installation
    and configuration instructions.

### Example

First, it is necessary to import `ase` and create the `ase.Atoms` object and attach
the desired calculator. This works exactly as in Python:

```@example ase
using PythonCall: pyimport, pylist

ase_build = pyimport("ase.build")
emt = pyimport("ase.calculators.emt")

h2 = ase_build.molecule("H2")
h2.calc = emt.EMT()
nothing # hide
```

Next, the [`ClassicalASEModel`](@ref) is created by passing the `ase.Atoms` object directly
to the model:

```@repl ase
using NQCModels
using NQCDInterfASE # This module contains all Python interface functionality. 
model = ClassicalASEModel(h2)
```

Now the model can be used in the same way as any of the previously introduced
analytic models.

```@repl ase
potential(model, rand(3, 2))
derivative(model, rand(3, 2))
```

!!! tip

    In theory, this should work with any of the ASE calculators that correctly implement
    the `get_potential_energy` and `get_forces` functions. For instance, you can use
    [SchNetPack (SPK)](https://github.com/atomistic-machine-learning/schnetpack) by
    passing their ASE calculator to the `ClassicalASEModel`.
    Take a look at [MLIP examples](@ref ml-pes-models) to learn more.

## MACE interface

[MACEModels.jl](https://github.com/NQCD/MACEModels.jl) is an interface to the [MACE](https://github.com/ACEsuit/mace) code. The package attempts to improve calculation speeds by directly interfacing with `MACE`, instead of going through `ase`.

More information on how to use MACE models is available in the [MACEModels.jl](@ref) API documentation. 

!!! tip

      To make the installation of a python environment to execute MACE in optional, the `MACEModel` is provided in a separate package, which needs to be installed with `]add MACEModels`. 

### Example

```julia
using NQCModels, NQCDynamics.jl, MACEModels

ensemble_paths = [
    "model01.model",
    "model02.model",
    "model03.model",
]

structure = Atoms([:N, :H, :H, :H]) # Define the structure
cell = InfiniteCell() # Set periodicity (none in this case)

# Create the MACE model
model = MACEModel(
    atoms, # Atoms (used for neighbour lists)
    cell, # Cell (used for neighbour lists)
    ensemble_paths; # Vector containing paths to one or more models
    device = "cpu", # Device (or vector of devices) to use
    default_dtype = Float32, # Data type of the models (Float32 or Float64)
    batch_size = 1, # Number of structures to evaluate at once (improves overall throughput)
)
```

## AtomsCalculators interface

[AtomsCalculators.jl](https://github.com/JuliaMolSim/AtomsCalculators.jl) provides a unified calculation interface for atomistic simulation engines within the [JuliaMolSim](https://github.com/JuliaMolSim) ecosystem.
NQCModels.jl supports bidirectional interoperability with AtomsCalculators.jl:

- **NQCModels → AtomsCalculators**: Any [`ClassicalModel`](@ref ClassicalModels.ClassicalModel) can be used directly as an AtomsCalculators calculator via `AtomsCalculators.potential_energy`, `AtomsCalculators.forces`, and `AtomsCalculators.virial`.
- **AtomsCalculators → NQCModels**: Any calculator implementing the AtomsCalculators interface can be wrapped with [`AtomsCalculatorsModel`](@ref) to produce a [`ClassicalModel`](@ref ClassicalModels.ClassicalModel) compatible with NQCDynamics.jl.

For detailed examples of both directions, see the [AtomsCalculators interoperability](@ref atomscalculators) page.

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
