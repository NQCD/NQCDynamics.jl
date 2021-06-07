# Overview

```@docs
NonadiabaticModels
```

## Abstract types

Models are categorised based upon the quantities that they provide.
Each of these abstract types defines a specific kind of model.
```@docs
NonadiabaticModels.AdiabaticModel
NonadiabaticModels.DiabaticModel
NonadiabaticModels.AdiabaticFrictionModel
NonadiabaticModels.DiabaticFrictionModel
```

## Implementing a new model
To add a new model the first step is to select the abstract type that it falls under.
After this, the new type should be created (see above for example)
and the following functions should be implemented:
```@docs
NonadiabaticModels.potential!
NonadiabaticModels.derivative!
NonadiabaticModels.friction!
```

!!! note

    It is required that the positions are provided to the model as a `DoFs*n_atoms` matrix.
    For a simple 1D model this is a `1*1` matrix,
    and for a 3D system this is a `3*N` matrix.

After the above non-allocating versions have been implemented,
allocating versions are automatically provided for testing and plotting convenience:
```@docs
NonadiabaticModels.potential
NonadiabaticModels.derivative
NonadiabaticModels.friction
```
