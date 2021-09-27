# Overview

`NonadiabaticModels.jl` is a package that provides a common interface for defining
models used for nonadiabatic dynamics simulations.
The interface is flexible enough such that concrete implementations can range from
simple 1D analytic forms to full-dimensional ab-initio models.

The advantage of this interface is that by implementing a common set of functions
it allows the same dynamics code to work for simple models and complex systems,
making scaling up as simple as possible.

The existing models in this package can be seen in the tree below.
Most of these are simple analytic models, with the complex specialist
models kept in separate packages.

```@example
import AbstractTrees
import InteractiveUtils: subtypes
import NonadiabaticModels: Model

AbstractTrees.children(x::Type) = subtypes(x)

AbstractTrees.print_tree(Model)
```

As is standard in Julia, only the types shown in the leaves of the tree are concrete
types, the branches are abstract types that group similar models together.
These abstract types denote the quantities provided by each model and the functions
that are implemented.

## Implementing a model

The simplest way to understand the `NonadiabticModels` interface is to walk through
an example implementation.

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
