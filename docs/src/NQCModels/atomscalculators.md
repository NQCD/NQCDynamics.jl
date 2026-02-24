```@setup logging
@info "Expanding src/NQCModels/atomscalculators.md..."
start_time = time()
```

# [AtomsCalculators interoperability](@id atomscalculators)

[AtomsCalculators.jl](https://github.com/JuliaMolSim/AtomsCalculators.jl) is a unified calculation interface for the [JuliaMolSim](https://github.com/JuliaMolSim) ecosystem.
NQCModels.jl supports bidirectional interoperability: any [`ClassicalModel`](@ref ClassicalModels.ClassicalModel) can be used as an AtomsCalculators calculator, and any AtomsCalculators calculator can be wrapped as a [`ClassicalModel`](@ref ClassicalModels.ClassicalModel).

## Using a `ClassicalModel` as an AtomsCalculators calculator

When AtomsCalculators.jl is loaded alongside NQCModels.jl, every [`ClassicalModel`](@ref ClassicalModels.ClassicalModel) automatically implements the AtomsCalculators interface.
This allows NQCModels to be used with any tool in the JuliaMolSim ecosystem that accepts an AtomsCalculators calculator, such as [Molly.jl](https://github.com/JuliaMolSim/Molly.jl) or geometry optimisers.

The three functions provided are:

- `AtomsCalculators.potential_energy(sys, model)` — returns the potential energy with Hartree units.
- `AtomsCalculators.forces(sys, model)` — returns atomic forces as a `Vector` of `SVector{3}` with Hartree/Bohr units.
- `AtomsCalculators.virial(sys, model)` — returns the virial stress tensor (zero for models that do not implement it).

```julia
using NQCModels
using AtomsCalculators
using AtomsBase
using Unitful, UnitfulAtomic

# Define a simple harmonic model
model = Harmonic()

# Build an AtomsBase system with one hydrogen atom
sys = isolated_system([
    :H => [0.0u"Å", 0.0u"Å", 0.5u"Å"],
])

# Evaluate energy and forces through the AtomsCalculators interface
E = AtomsCalculators.potential_energy(sys, model)
F = AtomsCalculators.forces(sys, model)
```

## Wrapping an AtomsCalculators calculator as a `ClassicalModel`

The [`AtomsCalculatorsModel`](@ref) type accepts any calculator that implements the AtomsCalculators interface and exposes it as a [`ClassicalModel`](@ref ClassicalModels.ClassicalModel), making it compatible with the full NQCDynamics.jl simulation framework.

The constructor requires the calculator object and either an `AtomsBase.AbstractSystem` or an `NQCBase.Structure` to define the atomic structure (species and cell):

```julia
using NQCModels
using NQCBase
using AtomsCalculators
using AtomsBase
using Unitful, UnitfulAtomic

# Assume `my_calc` is any AtomsCalculators-compatible calculator
# and `sys` is an AtomsBase.AbstractSystem describing the structure.

model = AtomsCalculatorsModel(my_calc, sys)

# The model can now be used like any other ClassicalModel
R = rand(3, length(sys))
E = potential(model, R)
D = derivative(model, R)
```

An `NQCBase.Structure` (containing `Atoms` and a `Cell`) can also be passed directly:

```julia
atoms = Atoms([:H, :H])
cell  = InfiniteCell()
structure = NQCBase.Structure(atoms, rand(3, 2), cell)

model = AtomsCalculatorsModel(my_calc, structure)
```

## Round-trip example

The bidirectional interface means that models can be converted in both directions while preserving numerical values.
The following example wraps a `Harmonic` model into an `AtomsCalculatorsModel` and verifies that the potential energies agree:

```julia
using NQCModels
using NQCBase
using AtomsCalculators

original_model = Harmonic()

atoms = Atoms([:H])
cell  = InfiniteCell()
R     = rand(3, 1)
structure = NQCBase.Structure(atoms, R, cell)

# Wrap the ClassicalModel as an AtomsCalculators calculator
# and then re-wrap it back as an AtomsCalculatorsModel.
roundtrip_model = AtomsCalculatorsModel(original_model, structure)

# The potential energies should agree to machine precision.
@assert isapprox(
    potential(original_model, R),
    potential(roundtrip_model, R),
    rtol = 1e-10,
)
```

```@setup logging
runtime = round(time() - start_time; digits=2)
@info "...done after $runtime s."
```
